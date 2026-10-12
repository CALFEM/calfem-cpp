//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// 2D beam elements, ported from calfem-python (core.py).

#include <calfem/elements/beam.h>

#include "detail/element_utils.h"

#include <cmath>
#include <stdexcept>

namespace calfem {

using namespace detail;

namespace {

// Local stiffness matrix of a 2D Euler-Bernoulli beam (axial + bending).

Matrix beam2LocalStiffness(double DEA, double DEI, double L)
{
    double L2 = L * L, L3 = L2 * L;

    Matrix K(6, 6);
    K << DEA / L, 0.0, 0.0, -DEA / L, 0.0, 0.0,                                         //
        0.0, 12.0 * DEI / L3, 6.0 * DEI / L2, 0.0, -12.0 * DEI / L3, 6.0 * DEI / L2,    //
        0.0, 6.0 * DEI / L2, 4.0 * DEI / L, 0.0, -6.0 * DEI / L2, 2.0 * DEI / L,        //
        -DEA / L, 0.0, 0.0, DEA / L, 0.0, 0.0,                                          //
        0.0, -12.0 * DEI / L3, -6.0 * DEI / L2, 0.0, 12.0 * DEI / L3, -6.0 * DEI / L2, //
        0.0, 6.0 * DEI / L2, 2.0 * DEI / L, 0.0, -6.0 * DEI / L2, 4.0 * DEI / L;
    return K;
}

// Local load vector for distributed loads qX (axial) and qY (transverse).

ColVec beam2LocalLoad(double qX, double qY, double L)
{
    ColVec f(6);
    f << qX / 2.0, qY / 2.0, qY * L / 12.0, qX / 2.0, qY / 2.0, -qY * L / 12.0;
    return L * f;
}

// Transforms local element matrices to global ones.

std::tuple<Matrix, ColVec> toGlobal(const Matrix &G, const Matrix &Kle, const ColVec &fle)
{
    return {G.transpose() * Kle * G, G.transpose() * fle};
}

const RowVec &zeroLoad2()
{
    static const RowVec zero = RowVec::Zero(2);
    return zero;
}

// Local displacements of a 2D beam split into axial (a1) and bending (a2) parts.

struct Beam2Local {
    Axis2D axis;
    ColVec edl;
    ColVec a2; // [v1, theta1, v2, theta2]
};

Beam2Local beam2Local(const RowVec &ex, const RowVec &ey, const RowVec &ed)
{
    Beam2Local b{axis2D(ex, ey), ColVec(), ColVec(4)};
    b.edl = beamTransform2D(b.axis) * ed.transpose();
    b.a2 << b.edl(1), b.edl(2), b.edl(4), b.edl(5);
    return b;
}

// Section forces [N, V, M] and displacements [u, v] for beam2s/beam2ws.

std::tuple<Matrix, Matrix, ColVec> beam2Section(const RowVec &ex, const RowVec &ey, double DEA, double DEI,
                                                const RowVec &ed, const RowVec &eq, double kX, double kY, int nep,
                                                const char *caller)
{
    auto local = beam2Local(ex, ey, ed);
    double L = local.axis.L;

    auto axial = axialSection(L, local.edl(0), local.edl(3), DEA, eq(0), kX, nep, caller);
    auto bending = bendingSection(L, local.a2, DEI, eq(1), kY, axial.X);

    Matrix es(nep, 3);
    es.col(0) = axial.N;
    es.col(1) = -DEI * bending.d3v;
    es.col(2) = DEI * bending.d2v;

    Matrix edi(nep, 2);
    edi.col(0) = axial.u;
    edi.col(1) = bending.v;

    return {es, edi, axial.X};
}

// Second order beam stiffness/load factors (beam2gxe).

struct SecondOrderFactors {
    double f2, f3, f4, f5, h;
};

SecondOrderFactors secondOrderFactors(double QX, double DEI, double L)
{
    const double eps = 1e-12;

    if (QX < -eps * DEI / (L * L))
    {
        double kL = std::sqrt(-QX / DEI) * L;
        double f1 = (kL / 2.0) / std::tan(kL / 2.0);
        double f2 = kL * kL / (12.0 * (1.0 - f1));
        double h = 6.0 * (2.0 / (kL * kL) - (1.0 + std::cos(kL)) / (kL * std::sin(kL)));
        return {f2, f1 / 4.0 + 3.0 * f2 / 4.0, -f1 / 2.0 + 3.0 * f2 / 2.0, f1 * f2, h};
    }
    else if (QX > eps * DEI / (L * L))
    {
        double kL = std::sqrt(QX / DEI) * L;
        double f1 = (kL / 2.0) / std::tanh(kL / 2.0);
        double f2 = -(1.0 / 12.0) * kL * kL / (1.0 - f1);
        double h = -6.0 * (2.0 / (kL * kL) - (1.0 + std::cosh(kL)) / (kL * std::sinh(kL)));
        return {f2, f1 / 4.0 + 3.0 * f2 / 4.0, -f1 / 2.0 + 3.0 * f2 / 2.0, f1 * f2, h};
    }

    return {1.0, 1.0, 1.0, 1.0, 1.0};
}

// Assembles es = [N, V, M] and edi = [u, v] for the second order beams, where the normal
// force includes the contribution of the transverse force through the slope.

std::tuple<Matrix, double, Matrix, ColVec> secondOrderResult(const AxialSection &axial, double DEI,
                                                             const Eigen::ArrayXd &v, const Eigen::ArrayXd &dv,
                                                             const Eigen::ArrayXd &d2v, const Eigen::ArrayXd &d3v)
{
    double QX = axial.N(0);

    Eigen::ArrayXd M = DEI * d2v;
    Eigen::ArrayXd V = -DEI * d3v;
    Eigen::ArrayXd N = QX + dv * V;

    Matrix es(axial.X.size(), 3);
    es.col(0) = N;
    es.col(1) = V;
    es.col(2) = M;

    Matrix edi(axial.X.size(), 2);
    edi.col(0) = axial.u;
    edi.col(1) = v;

    return {es, QX, edi, axial.X};
}

} // namespace

// ---- beam2e / beam2s ----

std::tuple<Matrix, ColVec> beam2e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &eq)
{
    double E = ep(0), A = ep(1), I = ep(2);
    auto axis = axis2D(ex, ey);

    return toGlobal(beamTransform2D(axis), beam2LocalStiffness(E * A, E * I, axis.L),
                    beam2LocalLoad(eq(0), eq(1), axis.L));
}

Matrix beam2e(const RowVec &ex, const RowVec &ey, const RowVec &ep)
{
    return std::get<0>(beam2e(ex, ey, ep, zeroLoad2()));
}

std::tuple<Matrix, Matrix, ColVec> beam2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                          const RowVec &eq, int nep)
{
    double E = ep(0), A = ep(1), I = ep(2);
    return beam2Section(ex, ey, E * A, E * I, ed, eq, 0.0, 0.0, nep, "beam2s");
}

Matrix beam2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, const RowVec &eq)
{
    return std::get<0>(beam2s(ex, ey, ep, ed, eq, 2));
}

Matrix beam2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed)
{
    return beam2s(ex, ey, ep, ed, zeroLoad2());
}

// ---- beam2we / beam2ws ----

std::tuple<Matrix, ColVec> beam2we(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &eq)
{
    double E = ep(0), A = ep(1), I = ep(2), kX = ep(3), kY = ep(4);
    auto axis = axis2D(ex, ey);
    double L = axis.L, L2 = L * L;

    Matrix Ks(6, 6);
    Ks << 140.0 * kX, 0.0, 0.0, 70.0 * kX, 0.0, 0.0,                          //
        0.0, 156.0 * kY, 22.0 * kY * L, 0.0, 54.0 * kY, -13.0 * kY * L,       //
        0.0, 22.0 * kY * L, 4.0 * kY * L2, 0.0, 13.0 * kY * L, -3.0 * kY * L2, //
        70.0 * kX, 0.0, 0.0, 140.0 * kX, 0.0, 0.0,                            //
        0.0, 54.0 * kY, 13.0 * kY * L, 0.0, 156.0 * kY, -22.0 * kY * L,       //
        0.0, -13.0 * kY * L, -3.0 * kY * L2, 0.0, -22.0 * kY * L, 4.0 * kY * L2;

    Matrix Kle = beam2LocalStiffness(E * A, E * I, L) + L / 420.0 * Ks;
    return toGlobal(beamTransform2D(axis), Kle, beam2LocalLoad(eq(0), eq(1), L));
}

Matrix beam2we(const RowVec &ex, const RowVec &ey, const RowVec &ep)
{
    return std::get<0>(beam2we(ex, ey, ep, zeroLoad2()));
}

std::tuple<Matrix, Matrix, ColVec> beam2ws(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                           const RowVec &eq, int nep)
{
    double E = ep(0), A = ep(1), I = ep(2), kX = ep(3), kY = ep(4);
    return beam2Section(ex, ey, E * A, E * I, ed, eq, kX, kY, nep, "beam2ws");
}

Matrix beam2ws(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, const RowVec &eq)
{
    return std::get<0>(beam2ws(ex, ey, ep, ed, eq, 2));
}

Matrix beam2ws(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed)
{
    return beam2ws(ex, ey, ep, ed, zeroLoad2());
}

// ---- beam2te / beam2ts (Timoshenko) ----

std::tuple<Matrix, ColVec> beam2te(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &eq)
{
    double E = ep(0), Gm = ep(1), A = ep(2), I = ep(3), ks = ep(4);
    double DEA = E * A, DEI = E * I, DGAks = Gm * A * ks;

    auto axis = axis2D(ex, ey);
    double L = axis.L, L2 = L * L, L3 = L2 * L;

    double m = (12.0 * DEI) / (L2 * DGAks);
    double f1 = 1.0 / (1.0 + m);
    double f2 = f1 * (1.0 + m / 4.0);
    double f3 = f1 * (1.0 - m / 2.0);

    Matrix Kle(6, 6);
    Kle << DEA / L, 0.0, 0.0, -DEA / L, 0.0, 0.0,                                                      //
        0.0, 12.0 * DEI * f1 / L3, 6.0 * DEI * f1 / L2, 0.0, -12.0 * DEI * f1 / L3, 6.0 * DEI * f1 / L2,    //
        0.0, 6.0 * DEI * f1 / L2, 4.0 * DEI * f2 / L, 0.0, -6.0 * DEI * f1 / L2, 2.0 * DEI * f3 / L,        //
        -DEA / L, 0.0, 0.0, DEA / L, 0.0, 0.0,                                                         //
        0.0, -12.0 * DEI * f1 / L3, -6.0 * DEI * f1 / L2, 0.0, 12.0 * DEI * f1 / L3, -6.0 * DEI * f1 / L2, //
        0.0, 6.0 * DEI * f1 / L2, 2.0 * DEI * f3 / L, 0.0, -6.0 * DEI * f1 / L2, 4.0 * DEI * f2 / L;

    return toGlobal(beamTransform2D(axis), Kle, beam2LocalLoad(eq(0), eq(1), L));
}

Matrix beam2te(const RowVec &ex, const RowVec &ey, const RowVec &ep)
{
    return std::get<0>(beam2te(ex, ey, ep, zeroLoad2()));
}

std::tuple<Matrix, Matrix, ColVec> beam2ts(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                           const RowVec &eq, int nep)
{
    double E = ep(0), Gm = ep(1), A = ep(2), I = ep(3), ks = ep(4);
    double DEA = E * A, DEI = E * I, DGAks = Gm * A * ks;
    double alpha = DEI / DGAks;
    double qX = eq(0), qY = eq(1);

    auto local = beam2Local(ex, ey, ed);
    double L = local.axis.L, L2 = L * L;

    auto axial = axialSection(L, local.edl(0), local.edl(3), DEA, qX, 0.0, nep, "beam2ts");

    Matrix C2(4, 4);
    C2 << L2 + 12.0 * alpha, 0.0, 0.0, 0.0,                                        //
        -12.0 * alpha / L, L2 + 6.0 * alpha, 12.0 * alpha / L, -6.0 * alpha,      //
        -3.0, -2.0 * L - 6.0 * alpha / L, 3.0, -L + 6.0 * alpha / L,              //
        2.0 / L, 1.0, -2.0 / L, 1.0;
    C2 *= 1.0 / (L2 + 12.0 * alpha);

    ColVec C2a = C2 * local.a2;
    double c0 = C2a(0), c1 = C2a(1), c2 = C2a(2), c3 = C2a(3);

    Eigen::ArrayXd x = axial.X.array();
    Eigen::ArrayXd x2 = x.square(), x3 = x2 * x;

    Eigen::ArrayXd v = c0 + c1 * x + c2 * x2 + c3 * x3;
    Eigen::ArrayXd dv = c1 + 2.0 * c2 * x + 3.0 * c3 * x2;
    Eigen::ArrayXd theta = c1 + 2.0 * c2 * x + (3.0 * x2 + 6.0 * alpha) * c3;
    Eigen::ArrayXd dtheta = 2.0 * c2 + 6.0 * c3 * x;

    if (DEI != 0.0)
    {
        v += (x2 * x2 - 2.0 * L * x3 + L2 * x2) * qY / (24.0 * DEI) + (-x2 / 2.0 + L * x / 2.0) * qY / DGAks;
        dv += (2.0 * x3 - 3.0 * L * x2 + L2 * x) * qY / (12.0 * DEI) + (-x + L / 2.0) * qY / DGAks;
        theta += (2.0 * x3 - 3.0 * L * x2 + L2 * x) * qY / (12.0 * DEI);
        dtheta += (6.0 * x2 - 6.0 * L * x + L2) * qY / (12.0 * DEI);
    }

    Matrix es(nep, 3);
    es.col(0) = axial.N;
    es.col(1) = DGAks * (dv - theta);
    es.col(2) = DEI * dtheta;

    Matrix edi(nep, 3);
    edi.col(0) = axial.u;
    edi.col(1) = v;
    edi.col(2) = theta;

    return {es, edi, axial.X};
}

Matrix beam2ts(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, const RowVec &eq)
{
    return std::get<0>(beam2ts(ex, ey, ep, ed, eq, 2));
}

Matrix beam2ts(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed)
{
    return beam2ts(ex, ey, ep, ed, zeroLoad2());
}

// ---- beam2ge / beam2gs (second order theory) ----

std::tuple<Matrix, ColVec> beam2ge(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX, double eq)
{
    double E = ep(0), A = ep(1), I = ep(2);
    auto axis = axis2D(ex, ey);
    double L = axis.L, L2 = L * L;
    double qY = eq;

    Matrix Ksle(6, 6);
    Ksle << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,         //
        0.0, 36.0, 3.0 * L, 0.0, -36.0, 3.0 * L,  //
        0.0, 3.0 * L, 4.0 * L2, 0.0, -3.0 * L, -L2, //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,             //
        0.0, -36.0, -3.0 * L, 0.0, 36.0, -3.0 * L, //
        0.0, 3.0 * L, -L2, 0.0, -3.0 * L, 4.0 * L2;

    Matrix Kle = beam2LocalStiffness(E * A, E * I, L) + QX / (30.0 * L) * Ksle;

    ColVec fle(6);
    fle << 0.0, 0.5, L / 12.0, 0.0, 0.5, -L / 12.0;
    fle *= qY * L;

    return toGlobal(beamTransform2D(axis), Kle, fle);
}

Matrix beam2ge(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX)
{
    return std::get<0>(beam2ge(ex, ey, ep, QX, 0.0));
}

std::tuple<Matrix, double, Matrix, ColVec> beam2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep,
                                                   const RowVec &ed, double QX, double eq, int nep)
{
    double E = ep(0), A = ep(1), I = ep(2);
    double DEA = E * A, DEI = E * I;

    auto local = beam2Local(ex, ey, ed);
    double L = local.axis.L, L2 = L * L, L3 = L2 * L;

    auto axial = axialSection(L, local.edl(0), local.edl(3), DEA, 0.0, 0.0, nep, "beam2gs");
    auto b = bendingSection(L, local.a2, DEI, eq, 0.0, axial.X);

    if (DEI != 0.0)
    {
        ColVec C2a = bendingInterpolation(L) * local.a2;
        double c2 = C2a(2), c3 = C2a(3);

        Eigen::ArrayXd x = axial.X.array();
        Eigen::ArrayXd x2 = x.square(), x3 = x2 * x, x4 = x3 * x, x5 = x4 * x;

        b.v += QX / DEI * ((x4 - 2.0 * L * x3 + L2 * x2) / 12.0 * c2 + (x5 - 3.0 * L2 * x3 + 2.0 * L3 * x2) / 20.0 * c3);
        b.dv += QX / DEI *
                ((2.0 * x3 - 3.0 * L * x2 + L2 * x) / 6.0 * c2 + (5.0 * x4 - 9.0 * L2 * x2 + 4.0 * L3 * x) / 20.0 * c3);
        b.d2v += QX / DEI * ((6.0 * x2 - 6.0 * L * x + L2) / 6.0 * c2 + (10.0 * x3 - 9.0 * L2 * x + 2.0 * L3) / 10.0 * c3);
        b.d3v += QX / DEI * ((2.0 * x - L) * c2 + (30.0 * x2 - 9.0 * L2) / 10.0 * c3);
    }

    return secondOrderResult(axial, DEI, b.v, b.dv, b.d2v, b.d3v);
}

std::tuple<Matrix, double> beam2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                   double QX, double eq)
{
    auto result = beam2gs(ex, ey, ep, ed, QX, eq, 2);
    return {std::get<0>(result), std::get<1>(result)};
}

std::tuple<Matrix, double> beam2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                   double QX)
{
    return beam2gs(ex, ey, ep, ed, QX, 0.0);
}

// ---- beam2gxe / beam2gxs (second order theory, exact solution) ----

std::tuple<Matrix, ColVec> beam2gxe(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX, double eq)
{
    double E = ep(0), A = ep(1), I = ep(2);
    double DEA = E * A, DEI = E * I;

    auto axis = axis2D(ex, ey);
    double L = axis.L, L2 = L * L, L3 = L2 * L;
    double qY = eq;

    auto f = secondOrderFactors(QX, DEI, L);

    Matrix Kle(6, 6);
    Kle << DEA / L, 0.0, 0.0, -DEA / L, 0.0, 0.0,                                                      //
        0.0, 12.0 * DEI * f.f5 / L3, 6.0 * DEI * f.f2 / L2, 0.0, -12.0 * DEI * f.f5 / L3, 6.0 * DEI * f.f2 / L2, //
        0.0, 6.0 * DEI * f.f2 / L2, 4.0 * DEI * f.f3 / L, 0.0, -6.0 * DEI * f.f2 / L2, 2.0 * DEI * f.f4 / L,     //
        -DEA / L, 0.0, 0.0, DEA / L, 0.0, 0.0,                                                         //
        0.0, -12.0 * DEI * f.f5 / L3, -6.0 * DEI * f.f2 / L2, 0.0, 12.0 * DEI * f.f5 / L3,             //
        -6.0 * DEI * f.f2 / L2,                                                                       //
        0.0, 6.0 * DEI * f.f2 / L2, 2.0 * DEI * f.f4 / L, 0.0, -6.0 * DEI * f.f2 / L2, 4.0 * DEI * f.f3 / L;

    ColVec fle(6);
    fle << 0.0, 0.5, L * f.h / 12.0, 0.0, 0.5, -L * f.h / 12.0;
    fle *= qY * L;

    return toGlobal(beamTransform2D(axis), Kle, fle);
}

Matrix beam2gxe(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX)
{
    return std::get<0>(beam2gxe(ex, ey, ep, QX, 0.0));
}

std::tuple<Matrix, double, Matrix, ColVec> beam2gxs(const RowVec &ex, const RowVec &ey, const RowVec &ep,
                                                    const RowVec &ed, double QX, double eq, int nep)
{
    double E = ep(0), A = ep(1), I = ep(2);
    double DEA = E * A, DEI = E * I;
    double qY = eq;

    auto local = beam2Local(ex, ey, ed);
    double L = local.axis.L;

    auto axial = axialSection(L, local.edl(0), local.edl(3), DEA, 0.0, 0.0, nep, "beam2gxs");

    Eigen::ArrayXd x = axial.X.array();
    auto n = x.size();
    Eigen::ArrayXd one = Eigen::ArrayXd::Ones(n);

    Eigen::ArrayXd v, dv, d2v, d3v;
    const double eps = 1e-12;

    if (QX < -eps * DEI / (L * L))
    {
        double k = std::sqrt(-QX / DEI);
        double kL = k * L;
        double s = std::sin(kL), c = std::cos(kL);

        Matrix C2(4, 4);
        C2 << k * (kL * s + c - 1.0), -kL * c + s, -k * (1.0 - c), -s + kL, //
            -(k * k) * s, -k * (1.0 - c), (k * k) * s, -k * (1.0 - c),      //
            -k * (1.0 - c), kL * c - s, k * (1.0 - c), s - kL,              //
            k * s, (kL * s + c - 1.0), -k * s, (1.0 - c);
        C2 *= 1.0 / (k * (-2.0 * (1.0 - c) + kL * s));

        ColVec C2a = C2 * local.a2;
        Eigen::ArrayXd ckx = (k * x).cos(), skx = (k * x).sin();

        v = C2a(0) + x * C2a(1) + ckx * C2a(2) + skx * C2a(3);
        dv = C2a(1) * one - k * skx * C2a(2) + k * ckx * C2a(3);
        d2v = -k * k * ckx * C2a(2) - k * k * skx * C2a(3);
        d3v = k * k * k * skx * C2a(2) - k * k * k * ckx * C2a(3);

        if (DEI != 0.0)
        {
            double r = (1.0 + c) / s;
            v += qY * std::pow(L, 4) / (2.0 * DEI) *
                 (r / (kL * kL * kL) * (-1.0 + ckx) + skx / (kL * kL * kL) + x * (-1.0 + x / L) / (kL * kL * L));
            dv += qY * std::pow(L, 3) / (2.0 * DEI) * (r / (kL * kL) * (-skx) + ckx / (kL * kL) + (-1.0 + 2.0 * x / L) / (kL * kL));
            d2v += qY * L * L / (2.0 * DEI) * (r / kL * (-ckx) - skx / kL + 2.0 / (kL * kL));
            d3v += qY * L / (2.0 * DEI) * (r * skx - ckx);
        }
    }
    else if (QX > eps * DEI / (L * L))
    {
        double k = std::sqrt(QX / DEI);
        double kL = k * L;
        double s = std::sinh(kL), c = std::cosh(kL);

        Matrix C2(4, 4);
        C2 << k * (-kL * s + c - 1.0), -kL * c + s, -k * (1.0 - c), -s + kL, //
            (k * k) * s, -k * (1.0 - c), -(k * k) * s, -k * (1.0 - c),       //
            -k * (1.0 - c), kL * c - s, k * (1.0 - c), s - kL,               //
            -k * s, (-kL * s + c - 1.0), k * s, (1.0 - c);
        C2 *= 1.0 / (k * (-2.0 * (1.0 - c) - kL * s));

        ColVec C2a = C2 * local.a2;
        Eigen::ArrayXd ckx = (k * x).cosh(), skx = (k * x).sinh();

        v = C2a(0) + x * C2a(1) + ckx * C2a(2) + skx * C2a(3);
        dv = C2a(1) * one + k * skx * C2a(2) + k * ckx * C2a(3);
        d2v = k * k * ckx * C2a(2) + k * k * skx * C2a(3);
        d3v = k * k * k * skx * C2a(2) + k * k * k * ckx * C2a(3);

        if (DEI != 0.0)
        {
            double r = (1.0 + c) / s;
            v += qY * std::pow(L, 4) / (2.0 * DEI) *
                 (r / (kL * kL * kL) * (-1.0 + ckx) - skx / (kL * kL * kL) + x * (1.0 - x / L) / (kL * kL * L));
            dv += qY * std::pow(L, 3) / (2.0 * DEI) * (r / (kL * kL) * skx - ckx / (kL * kL) + (1.0 - 2.0 * x / L) / (kL * kL));
            d2v += qY * L * L / (2.0 * DEI) * (r / kL * ckx - skx / kL - 2.0 / (kL * kL));
            d3v += qY * L / (2.0 * DEI) * (r * skx - ckx);
        }
    }
    else
    {
        auto b = bendingSection(L, local.a2, DEI, qY, 0.0, axial.X);
        v = b.v;
        dv = b.dv;
        d2v = b.d2v;
        d3v = b.d3v;
    }

    return secondOrderResult(axial, DEI, v, dv, d2v, d3v);
}

std::tuple<Matrix, double> beam2gxs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                    double QX, double eq)
{
    auto result = beam2gxs(ex, ey, ep, ed, QX, eq, 2);
    return {std::get<0>(result), std::get<1>(result)};
}

std::tuple<Matrix, double> beam2gxs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                    double QX)
{
    return beam2gxs(ex, ey, ep, ed, QX, 0.0);
}

// ---- beam2de / beam2ds (dynamics) ----

namespace {

struct Beam2Dynamic {
    Matrix Kle, Mle, Cle, G;
};

Beam2Dynamic beam2Dynamic(const RowVec &ex, const RowVec &ey, const RowVec &ep, const char *caller)
{
    if ((ep.size() != 4) && (ep.size() != 6))
        throw std::invalid_argument(std::string(caller) + ": ep must be [E, A, I, m] or [E, A, I, m, a, b].");

    double E = ep(0), A = ep(1), I = ep(2), m = ep(3);
    double a = ep.size() == 6 ? ep(4) : 0.0;
    double b = ep.size() == 6 ? ep(5) : 0.0;

    auto axis = axis2D(ex, ey);
    double L = axis.L, L2 = L * L;

    Beam2Dynamic d;
    d.Kle = beam2LocalStiffness(E * A, E * I, L);

    d.Mle.resize(6, 6);
    d.Mle << 140.0, 0.0, 0.0, 70.0, 0.0, 0.0,           //
        0.0, 156.0, 22.0 * L, 0.0, 54.0, -13.0 * L,     //
        0.0, 22.0 * L, 4.0 * L2, 0.0, 13.0 * L, -3.0 * L2, //
        70.0, 0.0, 0.0, 140.0, 0.0, 0.0,                //
        0.0, 54.0, 13.0 * L, 0.0, 156.0, -22.0 * L,     //
        0.0, -13.0 * L, -3.0 * L2, 0.0, -22.0 * L, 4.0 * L2;
    d.Mle *= m * L / 420.0;

    d.Cle = a * d.Mle + b * d.Kle;
    d.G = beamTransform2D(axis);
    return d;
}

} // namespace

std::tuple<Matrix, Matrix, Matrix> beam2de(const RowVec &ex, const RowVec &ey, const RowVec &ep)
{
    auto d = beam2Dynamic(ex, ey, ep, "beam2de");
    Matrix Gt = d.G.transpose();
    return {Gt * d.Kle * d.G, Gt * d.Mle * d.G, Gt * d.Cle * d.G};
}

Matrix beam2ds(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &ed, const Matrix &ev,
               const Matrix &ea)
{
    auto d = beam2Dynamic(ex, ey, ep, "beam2ds");

    Matrix es(ed.rows(), 6);

    for (Eigen::Index i = 0; i < ed.rows(); i++)
    {
        ColVec u = ed.row(i).transpose();
        ColVec v = ev.row(i).transpose();
        ColVec a = ea.row(i).transpose();
        es.row(i) = (d.Kle * d.G * u + d.Cle * d.G * v + d.Mle * d.G * a).transpose();
    }

    return es;
}

// ---- beam2crd ----

std::tuple<Matrix, Matrix> beam2crd(const Matrix &ex, const Matrix &ey, const Matrix &ed, double mag)
{
    const int nCoords = 21;
    auto nie = ed.rows();

    Matrix excd(nie, nCoords), eycd(nie, nCoords);

    for (Eigen::Index i = 0; i < nie; i++)
    {
        double bx = ex(i, 1) - ex(i, 0), by = ey(i, 1) - ey(i, 0);
        double L = std::sqrt(bx * bx + by * by);
        double n0 = bx / L, n1 = by / L;

        Axis2D axis{L, n0, n1};
        ColVec dl = beamTransform2D(axis) * ed.row(i).transpose();

        ColVec xl = ColVec::LinSpaced(nCoords, 0.0, L);
        ColVec one = ColVec::Ones(nCoords);

        // Axial displacement (linear)
        Matrix Cis(2, 2);
        Cis << -1.0, 1.0, //
            L, 0.0;
        Cis /= L;
        ColVec ds(2);
        ds << dl(0), dl(3);

        Matrix xlOne(nCoords, 2);
        xlOne << xl, one;
        ColVec ul = xlOne * Cis * ds;

        // Transverse displacement (cubic)
        Matrix Cib(4, 4);
        Cib << 12.0, 6.0 * L, -12.0, 6.0 * L,                    //
            -6.0 * L, -4.0 * L * L, 6.0 * L, -2.0 * L * L,       //
            0.0, L * L * L, 0.0, 0.0,                            //
            L * L * L, 0.0, 0.0, 0.0;
        Cib /= L * L * L;
        ColVec db(4);
        db << dl(1), dl(2), dl(4), dl(5);

        Matrix xlPoly(nCoords, 4);
        xlPoly << xl.array().cube().matrix() / 6.0, xl.array().square().matrix() / 2.0, xl, one;
        ColVec vl = xlPoly * Cib * db;

        // Back to global directions
        for (int k = 0; k < nCoords; k++)
        {
            double cx = n0 * ul(k) - n1 * vl(k);
            double cy = n1 * ul(k) + n0 * vl(k);
            excd(i, k) = n0 * xl(k) + ex(i, 0) + mag * cx;
            eycd(i, k) = n1 * xl(k) + ey(i, 0) + mag * cy;
        }
    }

    return {excd, eycd};
}

} // namespace calfem
