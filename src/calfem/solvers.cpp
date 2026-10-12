//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/solvers.h>

#include "detail/log.h"

#include <Eigen/Eigenvalues>

#include <string>
#include <vector>

namespace calfem {

namespace {

using detail::lastError;
using detail::logError;
using detail::logMessage;

[[noreturn]] void throwSolverError(const std::string &caller)
{
    throw SolverError(lastError().empty() ? caller + " failed." : lastError());
}

// Free/prescribed partitioning of the dofs of a system.

struct DofPartition {
    std::vector<Eigen::Index> ind;  // free dofs, zero-based, ascending
    std::vector<Eigen::Index> pind; // prescribed dofs, zero-based, in bcDofs order
};

bool checkSystem(const std::string &caller, Eigen::Index rows, Eigen::Index cols, Eigen::Index fSize)
{
    if (rows != cols)
    {
        logError(caller + ": K must be square (" + std::to_string(rows) + " x " + std::to_string(cols) + ").");
        return false;
    }

    if (fSize != rows)
    {
        logError(caller + ": f has " + std::to_string(fSize) + " rows, expected " + std::to_string(rows) + ".");
        return false;
    }

    return true;
}

bool partitionDofs(const std::string &caller, Eigen::Index nDofs, const IntColVec &bcDofs, const ColVec &bcVals,
                   DofPartition &p)
{
    if (bcDofs.size() != bcVals.size())
    {
        logError(caller + ": bcDofs and bcVals must have the same size (" + std::to_string(bcDofs.size()) + " and " +
            std::to_string(bcVals.size()) + ").");
        return false;
    }

    std::vector<bool> prescribed(static_cast<size_t>(nDofs), false);

    p.ind.clear();
    p.pind.clear();
    p.pind.reserve(static_cast<size_t>(bcDofs.size()));

    for (Eigen::Index i = 0; i < bcDofs.size(); i++)
    {
        Eigen::Index dof = bcDofs(i);

        if ((dof < 1) || (dof > nDofs))
        {
            logError(caller + ": prescribed dof " + std::to_string(dof) + " is outside 1.." + std::to_string(nDofs) + ".");
            return false;
        }

        if (prescribed[dof - 1])
        {
            logError(caller + ": dof " + std::to_string(dof) + " is prescribed more than once.");
            return false;
        }

        prescribed[dof - 1] = true;
        p.pind.push_back(dof - 1);
    }

    p.ind.reserve(static_cast<size_t>(nDofs) - p.pind.size());

    for (Eigen::Index i = 0; i < nDofs; i++)
        if (!prescribed[i])
            p.ind.push_back(i);

    return true;
}

// Splits K into Kff (free x free) and Kfp (free x prescribed). If Kp is given it receives
// the prescribed rows of K (prescribed x all).

void splitSparse(const SpMatrix &K, const DofPartition &p, SpMatrix &Kff, SpMatrix &Kfp, SpMatrix *Kp)
{
    auto nDofs = K.rows();
    auto nFree = static_cast<Eigen::Index>(p.ind.size());
    auto nPrescribed = static_cast<Eigen::Index>(p.pind.size());

    std::vector<Eigen::Index> freePos(static_cast<size_t>(nDofs), -1);
    std::vector<Eigen::Index> bcPos(static_cast<size_t>(nDofs), -1);

    for (Eigen::Index i = 0; i < nFree; i++)
        freePos[p.ind[i]] = i;

    for (Eigen::Index i = 0; i < nPrescribed; i++)
        bcPos[p.pind[i]] = i;

    TripletList ffList;
    TripletList fpList;
    TripletList pList;

    ffList.reserve(static_cast<size_t>(K.nonZeros()));

    for (Eigen::Index k = 0; k < K.outerSize(); ++k)
        for (SpMatrix::InnerIterator it(K, k); it; ++it)
        {
            auto r = freePos[it.row()];
            auto c = freePos[it.col()];

            if (r >= 0)
            {
                if (c >= 0)
                    ffList.emplace_back(r, c, it.value());
                else
                    fpList.emplace_back(r, bcPos[it.col()], it.value());
            }
            else if (Kp != nullptr)
                pList.emplace_back(bcPos[it.row()], it.col(), it.value());
        }

    Kff.resize(nFree, nFree);
    Kff.setFromTriplets(ffList.begin(), ffList.end());

    Kfp.resize(nFree, nPrescribed);
    Kfp.setFromTriplets(fpList.begin(), fpList.end());

    if (Kp != nullptr)
    {
        Kp->resize(nPrescribed, nDofs);
        Kp->setFromTriplets(pList.begin(), pList.end());
    }
}

// Warns (without failing) if the factorized system is close to singular, which usually means
// that the structure is not sufficiently supported. pivots is the diagonal of the Cholesky factor.

template <typename PivotVec, typename DiagVec>
void warnIfNearSingular(const std::string &caller, const PivotVec &pivots, const DiagVec &diagonal)
{
    if ((pivots.size() == 0) || (diagonal.size() == 0))
        return;

    double maxDiag = diagonal.cwiseAbs().maxCoeff();
    double minPivot = pivots.cwiseAbs().minCoeff();

    if (minPivot * minPivot <= 1e-12 * maxDiag)
        logMessage(caller + ": warning, K is nearly singular on the free dofs (relative pivot " +
                   std::to_string(minPivot * minPivot / maxDiag) + "). Check supports and mechanisms.");
}

std::string systemInfo(const std::string &caller, Eigen::Index nDofs, Eigen::Index nFree, Eigen::Index nonZeros)
{
    return caller + ": " + std::to_string(nDofs) + " dofs, " + std::to_string(nFree) + " free, " +
           std::to_string(nonZeros) + " non zeros in K.";
}

} // namespace

bool solveq(const Matrix &K, const ColVec &f, const IntColVec &bcDofs, const ColVec &bcVals, ColVec &a, ColVec &Q)
{
    if (!checkSystem("solveq", K.rows(), K.cols(), f.size()))
        return false;

    auto nDofs = K.rows();

    DofPartition p;

    if (!partitionDofs("solveq", nDofs, bcDofs, bcVals, p))
        return false;

    a.resize(nDofs);
    a.setZero();

    if (!p.ind.empty())
    {
        Matrix Ksys = K(p.ind, p.ind);
        ColVec fsys = f(p.ind) - K(p.ind, p.pind) * bcVals; // [ ind x 1 ] - [ ind x pind ] x [ pind ]

        Eigen::LLT<Matrix> llt(Ksys);

        if (llt.info() != Eigen::Success)
        {
            logError("solveq: factorization failed, K is not positive definite on the free dofs.");
            return false;
        }

        warnIfNearSingular("solveq", llt.matrixLLT().diagonal(), Ksys.diagonal());

        ColVec asys = llt.solve(fsys);
        a(p.ind) = asys;
    }

    a(p.pind) = bcVals;

    Q = K * a - f;

    return true;
}

bool spsolveq(const SpMatrix &K, const ColVec &f, const IntColVec &bcDofs, const ColVec &bcVals, ColVec &a, ColVec &Q)
{
    if (!checkSystem("spsolveq", K.rows(), K.cols(), f.size()))
        return false;

    auto nDofs = K.rows();

    DofPartition p;

    if (!partitionDofs("spsolveq", nDofs, bcDofs, bcVals, p))
        return false;

    SpMatrix Ksys;
    SpMatrix Ksysf;

    splitSparse(K, p, Ksys, Ksysf, nullptr);

    logMessage(systemInfo("spsolveq", nDofs, Ksys.rows(), K.nonZeros()));

    a.resize(nDofs);
    a.setZero();

    if (!p.ind.empty())
    {
        Eigen::SimplicialLLT<SpMatrix> solver;
        solver.compute(Ksys);

        if (solver.info() != Eigen::Success)
        {
            logError("spsolveq: factorization failed, K is not positive definite on the free dofs.");
            return false;
        }

        warnIfNearSingular("spsolveq", solver.matrixL().nestedExpression().diagonal(), Ksys.diagonal());

        ColVec fsys = f(p.ind) - Ksysf * bcVals;
        ColVec asys = solver.solve(fsys);

        if (solver.info() != Eigen::Success)
        {
            logError("spsolveq: solve failed.");
            return false;
        }

        a(p.ind) = asys;
    }

    a(p.pind) = bcVals;

    Q = K * a - f;

    return true;
}

// ---- SparseSolver ----

SparseSolver::SparseSolver()
{}

std::shared_ptr<SparseSolver> SparseSolver::create()
{
    return std::make_shared<SparseSolver>();
}

bool SparseSolver::setup(const SpMatrix &K, const IntColVec &bcDofs, const ColVec &bcVals)
{
    m_isSetup = false;
    m_factorized = false;

    if (!checkSystem("SparseSolver::setup", K.rows(), K.cols(), K.rows()))
        return false;

    DofPartition p;

    if (!partitionDofs("SparseSolver::setup", K.rows(), bcDofs, bcVals, p))
        return false;

    splitSparse(K, p, m_Ksys, m_Ksysf, &m_Kp);

    m_ind = std::move(p.ind);
    m_pind = std::move(p.pind);
    m_bcVals = bcVals;

    m_nDofs = K.rows();
    m_nFreeDofs = static_cast<Eigen::Index>(m_ind.size());

    logMessage(systemInfo("SparseSolver::setup", m_nDofs, m_nFreeDofs, K.nonZeros()));

    m_isSetup = true;

    return true;
}

bool SparseSolver::factorize()
{
    m_factorized = false;

    if (m_nFreeDofs > 0)
    {
        m_solver.compute(m_Ksys);

        if (m_solver.info() != Eigen::Success)
        {
            logError("SparseSolver: factorization failed, K is not positive definite on the free dofs.");
            return false;
        }

        warnIfNearSingular("SparseSolver", m_solver.matrixL().nestedExpression().diagonal(), m_Ksys.diagonal());
    }

    m_factorized = true;

    return true;
}

bool SparseSolver::solveSystem(const ColVec &f, ColVec &a, ColVec &Q)
{
    if (f.size() != m_nDofs)
    {
        logError("SparseSolver: f has " + std::to_string(f.size()) + " rows, expected " + std::to_string(m_nDofs) + ".");
        return false;
    }

    a.resize(m_nDofs);
    a.setZero();

    Q.resize(m_nDofs);
    Q.setZero();

    if (m_nFreeDofs > 0)
    {
        ColVec fsys = f(m_ind) - m_Ksysf * m_bcVals;
        ColVec asys = m_solver.solve(fsys);

        if (m_solver.info() != Eigen::Success)
        {
            logError("SparseSolver: solve failed.");
            return false;
        }

        a(m_ind) = asys;

        // Q = K a - f on the free dofs (equilibrium residual)

        ColVec Qf = m_Ksys * asys + m_Ksysf * m_bcVals - f(m_ind);
        Q(m_ind) = Qf;
    }

    a(m_pind) = m_bcVals;

    // Q = K a - f on the prescribed dofs (reaction forces)

    if (!m_pind.empty())
    {
        ColVec Qp = m_Kp * a - f(m_pind);
        Q(m_pind) = Qp;
    }

    return true;
}

bool SparseSolver::solve(const ColVec &f, ColVec &a, ColVec &Q)
{
    if (!m_isSetup)
    {
        logError("SparseSolver::solve: setup() has not completed successfully.");
        return false;
    }

    if (!factorize())
        return false;

    return solveSystem(f, a, Q);
}

bool SparseSolver::recompute(const ColVec &f, ColVec &a, ColVec &Q)
{
    if (!m_isSetup)
    {
        logError("SparseSolver::recompute: setup() has not completed successfully.");
        return false;
    }

    if (!m_factorized && !factorize())
        return false;

    return solveSystem(f, a, Q);
}

bool SparseSolver::isFactorized() const
{
    return m_factorized;
}

Eigen::Index SparseSolver::nFreeDofs() const
{
    return m_nFreeDofs;
}

// ---- Value-returning versions (calfem-python style) ----

std::tuple<ColVec, ColVec> solveq(const Matrix &K, const ColVec &f, const IntColVec &bcDofs)
{
    return solveq(K, f, bcDofs, ColVec::Zero(bcDofs.size()));
}

std::tuple<ColVec, ColVec> solveq(const Matrix &K, const ColVec &f, const IntColVec &bcDofs, const ColVec &bcVals)
{
    ColVec a, r;
    lastError().clear();

    if (!solveq(K, f, bcDofs, bcVals, a, r))
        throwSolverError("solveq");

    return {a, r};
}

std::tuple<ColVec, ColVec> spsolveq(const SpMatrix &K, const ColVec &f, const IntColVec &bcDofs)
{
    return spsolveq(K, f, bcDofs, ColVec::Zero(bcDofs.size()));
}

std::tuple<ColVec, ColVec> spsolveq(const SpMatrix &K, const ColVec &f, const IntColVec &bcDofs,
                                    const ColVec &bcVals)
{
    ColVec a, r;
    lastError().clear();

    if (!spsolveq(K, f, bcDofs, bcVals, a, r))
        throwSolverError("spsolveq");

    return {a, r};
}

bool SparseSolver::setup(const SpMatrix &K, const IntColVec &bcDofs)
{
    return setup(K, bcDofs, ColVec::Zero(bcDofs.size()));
}

std::tuple<ColVec, ColVec> SparseSolver::solve(const ColVec &f)
{
    ColVec a, r;
    lastError().clear();

    if (!solve(f, a, r))
        throwSolverError("SparseSolver::solve");

    return {a, r};
}

std::tuple<ColVec, ColVec> SparseSolver::recompute(const ColVec &f)
{
    ColVec a, r;
    lastError().clear();

    if (!recompute(f, a, r))
        throwSolverError("SparseSolver::recompute");

    return {a, r};
}

// ---- Reduction, static condensation and eigenvalues ----

namespace {

// Zero-based dofs of an ndof system that are not in the 1-based list b (ascending).
std::vector<Eigen::Index> remainingDofs(const std::string &caller, Eigen::Index nDofs, const IntColVec &b)
{
    std::vector<bool> removed(static_cast<size_t>(nDofs), false);

    for (Eigen::Index i = 0; i < b.size(); i++)
    {
        if ((b(i) < 1) || (b(i) > nDofs))
        {
            logError(caller + ": dof " + std::to_string(b(i)) + " is outside 1.." + std::to_string(nDofs) + ".");
            throwSolverError(caller);
        }
        removed[b(i) - 1] = true;
    }

    std::vector<Eigen::Index> remaining;
    for (Eigen::Index i = 0; i < nDofs; i++)
        if (!removed[i])
            remaining.push_back(i);

    return remaining;
}

std::vector<Eigen::Index> zeroBased(const IntColVec &dofs)
{
    std::vector<Eigen::Index> idx;
    for (Eigen::Index i = 0; i < dofs.size(); i++)
        idx.push_back(dofs(i) - 1);
    return idx;
}

} // namespace

Matrix red(const Matrix &A, const IntColVec &b)
{
    auto fdof = remainingDofs("red", A.rows(), b);
    return A(fdof, fdof);
}

Matrix red(const Matrix &A, const Matrix &bc)
{
    IntColVec b(bc.rows());
    for (Eigen::Index i = 0; i < bc.rows(); i++)
        b(i) = static_cast<int>(bc(i, 0));
    return red(A, b);
}

std::tuple<Matrix, ColVec> statcon(const Matrix &K, const ColVec &f, const IntColVec &cd)
{
    auto aindx = remainingDofs("statcon", K.rows(), cd);
    auto bindx = zeroBased(cd);

    Matrix Kaa = K(aindx, aindx);
    Matrix Kab = K(aindx, bindx);
    Matrix Kbb = K(bindx, bindx);

    ColVec fa = f(aindx);
    ColVec fb = f(bindx);

    Eigen::FullPivLU<Matrix> lu(Kbb);

    if (!lu.isInvertible())
    {
        logError("statcon: the condensed part of K is singular.");
        throwSolverError("statcon");
    }

    Matrix KbbInv = lu.inverse();

    Matrix K1 = Kaa - Kab * KbbInv * Kab.transpose();
    ColVec f1 = fa - Kab * KbbInv * fb;

    return {K1, f1};
}

std::tuple<ColVec, Matrix> eigen(const Matrix &K, const Matrix &M, const IntColVec &b)
{
    if ((K.rows() != K.cols()) || (M.rows() != K.rows()) || (M.cols() != K.cols()))
    {
        logError("eigen: K and M must be square and of the same size.");
        throwSolverError("eigen");
    }

    auto fdof = remainingDofs("eigen", K.rows(), b);

    Matrix Kf = K(fdof, fdof);
    Matrix Mf = M(fdof, fdof);

    // GeneralizedSelfAdjointEigenSolver does not report a failed Cholesky factorization of M,
    // so check positive definiteness first.
    if (Eigen::LLT<Matrix>(Mf).info() != Eigen::Success)
    {
        logError("eigen: M must be symmetric positive definite on the free dofs.");
        throwSolverError("eigen");
    }

    // Eigenvalues ascending, eigenvectors normalized so that X^T M X = I
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(Kf, Mf);

    if (solver.info() != Eigen::Success)
    {
        logError("eigen: the eigenvalue problem could not be solved.");
        throwSolverError("eigen");
    }

    Matrix X = Matrix::Zero(K.rows(), static_cast<Eigen::Index>(fdof.size()));
    X(fdof, Eigen::all) = solver.eigenvectors();

    return {solver.eigenvalues(), X};
}

std::tuple<ColVec, Matrix> eigen(const Matrix &K, const Matrix &M)
{
    return eigen(K, M, IntColVec(0));
}

} // namespace calfem
