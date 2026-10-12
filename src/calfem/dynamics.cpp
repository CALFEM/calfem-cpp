//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Time integration (step1, step2) and time functions (gfunc), ported from calfem-python.

#include <calfem/dynamics.h>
#include <calfem/solvers.h>

#include "detail/log.h"

#include <Eigen/LU>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace calfem {

namespace {

using detail::logError;

[[noreturn]] void fail(const std::string &message)
{
    logError(message);
    throw SolverError(message);
}

// calfem-python decides whether f, bc and times are given with np.array(x).any().
bool anyNonZero(const Matrix &m)
{
    return (m.size() > 0) && (m.array() != 0.0).any();
}

// Factorization of a dense system matrix that must be non-singular.
Eigen::PartialPivLU<Matrix> factorize(const Matrix &A, const std::string &message)
{
    Eigen::PartialPivLU<Matrix> lu(A);

    if (!(lu.rcond() > std::numeric_limits<double>::epsilon()))
        fail(message);

    return lu;
}

// Time stepping setup shared by step1 and step2.

struct StepSetup {
    Eigen::Index ndof = 0;
    int nstep = 0;
    Matrix tf; // load, one column per step (0..nstep)

    bool bound = false;
    Matrix pa;  // prescribed values, one column per step
    Matrix pda; // prescribed first derivatives
    std::vector<Eigen::Index> pdof; // zero-based, in bc order
    std::vector<Eigen::Index> fdof; // zero-based, ascending

    bool storeAll = true; // store every step in modelhist (no times given)
    std::vector<double> times;
    std::vector<Eigen::Index> dofs; // zero-based dofs for dofhist
};

StepSetup stepSetup(const char *caller, const Matrix &K, const Matrix &f, const ColVec &a0, const Matrix &bc,
                    double dt, double tottime, const RowVec &times, const IntColVec &dofs)
{
    std::string name(caller);
    StepSetup s;
    s.ndof = K.rows();

    if ((K.cols() != s.ndof) || (a0.size() != s.ndof))
        fail(name + ": K must be square and a0 must have one value per dof.");

    if (!(dt > 0.0))
        fail(name + ": the time step dt must be positive.");

    // Number of steps: from the columns of f or bc if they give a time history, limited by T/dt
    int nstep = 1;

    bool hasLoad = anyNonZero(f);
    if (hasLoad)
    {
        if (f.rows() != s.ndof)
            fail(name + ": f must have one row per dof.");
        if (f.cols() > 1)
            nstep = static_cast<int>(f.cols()) - 1;
    }

    s.bound = anyNonZero(bc);
    if (s.bound)
    {
        if (bc.cols() < 2)
            fail(name + ": bc must have at least two columns [dof, value].");
        if (bc.cols() > 2)
            nstep = static_cast<int>(bc.cols()) - 2;
    }

    int ns = static_cast<int>(tottime / dt);
    if ((ns < nstep) || (nstep == 1))
        nstep = ns;

    s.nstep = nstep;

    // Load history
    s.tf = Matrix::Zero(s.ndof, nstep + 1);
    if (hasLoad)
    {
        if (f.cols() == 1)
            s.tf = f.col(0).replicate(1, nstep + 1);
        else if (f.cols() >= nstep + 1)
            s.tf = f;
        else
            fail(name + ": f has fewer time steps than bc.");
    }

    // Prescribed values
    if (s.bound)
    {
        auto nrb = bc.rows();

        if (bc.cols() == 2)
        {
            s.pa = bc.col(1).replicate(1, nstep + 1);
            s.pda = Matrix::Zero(nrb, nstep + 1);
        }
        else
        {
            s.pa = bc.rightCols(bc.cols() - 1);
            s.pda.resize(nrb, s.pa.cols());
            s.pda.col(0) = (s.pa.col(1) - s.pa.col(0)) / dt;
            s.pda.rightCols(s.pa.cols() - 1) =
                (s.pa.rightCols(s.pa.cols() - 1) - s.pa.leftCols(s.pa.cols() - 1)) / dt;
        }

        std::vector<bool> prescribed(static_cast<size_t>(s.ndof), false);

        for (Eigen::Index i = 0; i < nrb; i++)
        {
            auto dof = static_cast<Eigen::Index>(bc(i, 0));
            if ((dof < 1) || (dof > s.ndof))
                fail(name + ": prescribed dof " + std::to_string(dof) + " is outside 1.." + std::to_string(s.ndof) +
                     ".");
            s.pdof.push_back(dof - 1);
            prescribed[dof - 1] = true;
        }

        for (Eigen::Index i = 0; i < s.ndof; i++)
            if (!prescribed[i])
                s.fdof.push_back(i);
    }
    else
    {
        for (Eigen::Index i = 0; i < s.ndof; i++)
            s.fdof.push_back(i);
    }

    // Storage
    s.storeAll = !anyNonZero(times);
    if (!s.storeAll)
        s.times.assign(times.data(), times.data() + times.size());

    // calfem-python uses dofhist only if all dofs are non-zero (np.array(dofs).all())
    if ((dofs.size() > 0) && (dofs.array() != 0).all())
    {
        for (Eigen::Index i = 0; i < dofs.size(); i++)
        {
            if ((dofs(i) < 1) || (dofs(i) > s.ndof))
                fail(name + ": history dof " + std::to_string(dofs(i)) + " is outside 1.." + std::to_string(s.ndof) +
                     ".");
            s.dofs.push_back(dofs(i) - 1);
        }
    }

    return s;
}

// Records results in modelhist and dofhist.

class HistoryRecorder {
public:
    HistoryRecorder(const StepSetup &s, bool secondOrder) : m_s(s), m_secondOrder(secondOrder)
    {
        auto nCols = s.storeAll ? s.nstep + 1 : static_cast<Eigen::Index>(s.times.size());
        auto nDofs = static_cast<Eigen::Index>(s.dofs.size());

        m_model.a = Matrix::Zero(s.ndof, nCols);
        m_model.da = Matrix::Zero(s.ndof, nCols);

        if (nDofs > 0)
        {
            m_dof.a = Matrix::Zero(nDofs, s.nstep + 1);
            m_dof.da = Matrix::Zero(nDofs, s.nstep + 1);
        }

        if (secondOrder)
        {
            m_model.d2a = Matrix::Zero(s.ndof, nCols);
            if (nDofs > 0)
                m_dof.d2a = Matrix::Zero(nDofs, s.nstep + 1);
        }
    }

    // Values at t = 0 (stored for the times option only if the first time is 0, as in calfem-python)
    void initial(const ColVec &a, const ColVec &da, const ColVec &d2a)
    {
        if (m_s.storeAll)
            store(m_model, 0, a, da, d2a);
        else if (!m_s.times.empty() && (m_s.times[0] == 0.0))
            store(m_model, m_itime++, a, da, d2a);

        storeDofs(0, a, da, d2a);
    }

    void step(int j, double time, const ColVec &a, const ColVec &da, const ColVec &d2a)
    {
        if (m_s.storeAll)
            store(m_model, j, a, da, d2a);
        else if ((m_itime < m_s.times.size()) && (time >= m_s.times[m_itime]))
            store(m_model, m_itime++, a, da, d2a);

        storeDofs(j, a, da, d2a);
    }

    std::tuple<TimeHistory, TimeHistory> result() const
    {
        return {m_model, m_dof};
    }

private:
    void store(TimeHistory &h, Eigen::Index col, const ColVec &a, const ColVec &da, const ColVec &d2a)
    {
        h.a.col(col) = a;
        h.da.col(col) = da;
        if (m_secondOrder)
            h.d2a.col(col) = d2a;
    }

    void storeDofs(Eigen::Index col, const ColVec &a, const ColVec &da, const ColVec &d2a)
    {
        if (m_s.dofs.empty())
            return;

        m_dof.a.col(col) = a(m_s.dofs);
        m_dof.da.col(col) = da(m_s.dofs);
        if (m_secondOrder)
            m_dof.d2a.col(col) = d2a(m_s.dofs);
    }

    const StepSetup &m_s;
    bool m_secondOrder;
    size_t m_itime = 0;
    TimeHistory m_model;
    TimeHistory m_dof;
};

} // namespace

// ---- step1 ----

std::tuple<TimeHistory, TimeHistory> step1(const Matrix &K, const Matrix &C, const Matrix &f, const ColVec &a0,
                                           const Matrix &bc, const RowVec &ip, const RowVec &times,
                                           const IntColVec &dofs)
{
    if (ip.size() != 3)
        fail("step1: ip must be [dt, T, alpha].");

    double dt = ip(0), tottime = ip(1), alpha = ip(2);
    double a1 = (1.0 - alpha) * dt;
    double a2 = alpha * dt;

    auto s = stepSetup("step1", K, f, a0, bc, dt, tottime, times, dofs);

    if ((C.rows() != s.ndof) || (C.cols() != s.ndof))
        fail("step1: C must have the same size as K.");

    auto Clu = factorize(C, "step1: C is singular.");
    ColVec da0 = Clu.solve(ColVec(s.tf.col(0) - K * a0));

    HistoryRecorder recorder(s, false);
    recorder.initial(a0, da0, ColVec());

    const auto &fd = s.fdof;
    const auto &pd = s.pdof;

    Matrix Kff = s.bound ? Matrix(K(fd, fd)) : K;
    Matrix Kfp = s.bound ? Matrix(K(fd, pd)) : Matrix();
    Matrix Cfp = s.bound ? Matrix(C(fd, pd)) : Matrix();
    Matrix Keff = s.bound ? Matrix(C(fd, fd) + a2 * Kff) : Matrix(C + a2 * K);

    auto lu = factorize(Keff, "step1: the effective system matrix is singular.");

    ColVec anew = a0(fd);
    ColVec danew = da0(fd);

    ColVec tempa = ColVec::Zero(s.ndof);
    ColVec tempda = ColVec::Zero(s.ndof);

    for (int j = 1; j <= s.nstep; j++)
    {
        double time = dt * j;

        ColVec apred = anew + a1 * danew;
        ColVec reff;

        if (s.bound)
        {
            ColVec pdeff = Cfp * s.pda.col(j) + Kfp * s.pa.col(j);
            reff = s.tf.col(j)(fd) - Kff * apred - pdeff;
        }
        else
            reff = s.tf.col(j) - K * apred;

        danew = lu.solve(reff);
        anew = apred + a2 * danew;

        if (s.bound)
        {
            tempa(pd) = s.pa.col(j);
            tempda(pd) = s.pda.col(j);
        }

        tempa(fd) = anew;
        tempda(fd) = danew;

        recorder.step(j, time, tempa, tempda, ColVec());
    }

    return recorder.result();
}

// ---- step2 ----

std::tuple<TimeHistory, TimeHistory> step2(const Matrix &K, const Matrix &C, const Matrix &M, const Matrix &f,
                                           const ColVec &a0, const ColVec &da0, const Matrix &bc, const RowVec &ip,
                                           const RowVec &times, const IntColVec &dofs)
{
    if (ip.size() != 4)
        fail("step2: ip must be [dt, T, alpha, delta].");

    double dt = ip(0), tottime = ip(1), alpha = ip(2), delta = ip(3);
    double b1 = dt * dt * 0.5 * (1.0 - 2.0 * alpha);
    double b2 = (1.0 - delta) * dt;
    double b3 = delta * dt;
    double b4 = alpha * dt * dt;

    auto s = stepSetup("step2", K, f, a0, bc, dt, tottime, times, dofs);

    if ((M.rows() != s.ndof) || (M.cols() != s.ndof))
        fail("step2: M must have the same size as K.");

    if (da0.size() != s.ndof)
        fail("step2: da0 must have one value per dof.");

    // No damping when C is empty or zero
    Matrix Cd = anyNonZero(C) ? C : Matrix(Matrix::Zero(s.ndof, s.ndof));

    if ((Cd.rows() != s.ndof) || (Cd.cols() != s.ndof))
        fail("step2: C must have the same size as K.");

    auto Mlu = factorize(M, "step2: M is singular.");
    ColVec d2a0 = Mlu.solve(ColVec(s.tf.col(0) - Cd * da0 - K * a0));

    HistoryRecorder recorder(s, true);
    recorder.initial(a0, da0, d2a0);

    const auto &fd = s.fdof;
    const auto &pd = s.pdof;

    Matrix Kff = s.bound ? Matrix(K(fd, fd)) : K;
    Matrix Cff = s.bound ? Matrix(Cd(fd, fd)) : Cd;
    Matrix Kfp = s.bound ? Matrix(K(fd, pd)) : Matrix();
    Matrix Cfp = s.bound ? Matrix(Cd(fd, pd)) : Matrix();
    Matrix Keff = s.bound ? Matrix(M(fd, fd) + b3 * Cff + b4 * Kff) : Matrix(M + b3 * Cd + b4 * K);

    auto lu = factorize(Keff, "step2: the effective system matrix is singular.");

    ColVec anew = a0(fd);
    ColVec danew = da0(fd);
    ColVec d2anew = d2a0(fd);

    ColVec tempa = ColVec::Zero(s.ndof);
    ColVec tempda = ColVec::Zero(s.ndof);
    ColVec tempd2a = ColVec::Zero(s.ndof);

    for (int j = 1; j <= s.nstep; j++)
    {
        double time = dt * j;

        ColVec apred = anew + dt * danew + b1 * d2anew;
        ColVec dapred = danew + b2 * d2anew;
        ColVec reff;

        if (s.bound)
        {
            ColVec pdeff = Cfp * s.pda.col(j) + Kfp * s.pa.col(j);
            reff = s.tf.col(j)(fd) - Cff * dapred - Kff * apred - pdeff;
        }
        else
            reff = s.tf.col(j) - Cd * dapred - K * apred;

        d2anew = lu.solve(reff);
        anew = apred + b4 * d2anew;
        danew = dapred + b3 * d2anew;

        if (s.bound)
        {
            tempa(pd) = s.pa.col(j);
            tempda(pd) = s.pda.col(j);
        }

        tempa(fd) = anew;
        tempda(fd) = danew;
        tempd2a(fd) = d2anew;

        recorder.step(j, time, tempa, tempda, tempd2a);
    }

    return recorder.result();
}

// ---- gfunc ----

std::tuple<ColVec, ColVec> gfunc(const Matrix &G, double dt)
{
    if ((G.rows() < 1) || (G.cols() < 2))
        fail("gfunc: G must have at least one row [t, g].");

    if (!(dt > 0.0))
        fail("gfunc: the time step dt must be positive.");

    double t0 = G(0, 0);
    double tend = G(G.rows() - 1, 0);

    auto n = static_cast<Eigen::Index>(std::floor((tend - t0) / dt + 1e-9)) + 1;

    ColVec ti(n), g(n);

    for (Eigen::Index i = 0; i < n; i++)
    {
        double t = t0 + i * dt;
        ti(i) = t;

        // Linear interpolation, constant outside [t0, tend] (as np.interp)
        if (t <= G(0, 0))
            g(i) = G(0, 1);
        else if (t >= tend)
            g(i) = G(G.rows() - 1, 1);
        else
        {
            Eigen::Index k = 0;
            while (G(k + 1, 0) <= t)
                k++;

            double slope = (G(k + 1, 1) - G(k, 1)) / (G(k + 1, 0) - G(k, 0));
            g(i) = slope * (t - G(k, 0)) + G(k, 1);
        }
    }

    return {ti, g};
}

} // namespace calfem
