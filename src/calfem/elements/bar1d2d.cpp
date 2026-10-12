//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// 1D and 2D bar elements, ported from calfem-python (core.py).

#include <calfem/elements/bar.h>

#include "detail/element_utils.h"

#include <cmath>

namespace calfem {

using namespace detail;

namespace {

double length1D(const RowVec &ex)
{
    return std::abs(ex(1) - ex(0));
}

Matrix axialStiffness(double DEA, double L)
{
    Matrix K(2, 2);
    K << 1.0, -1.0, //
        -1.0, 1.0;
    return DEA / L * K;
}

ColVec axialLoad(double qX, double L)
{
    ColVec f(2);
    f << 0.5, 0.5;
    return qX * L * f;
}

} // namespace

// ---- bar1e / bar1s ----

std::tuple<Matrix, ColVec> bar1e(const RowVec &ex, const RowVec &ep, double eq)
{
    double DEA = ep(0) * ep(1);
    double L = length1D(ex);
    return {axialStiffness(DEA, L), axialLoad(eq, L)};
}

Matrix bar1e(const RowVec &ex, const RowVec &ep)
{
    return std::get<0>(bar1e(ex, ep, 0.0));
}

std::tuple<ColVec, ColVec, ColVec> bar1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep)
{
    double DEA = ep(0) * ep(1);
    auto s = axialSection(length1D(ex), ed(0), ed(1), DEA, eq, 0.0, nep, "bar1s");
    return {s.N, s.u, s.X};
}

ColVec bar1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq)
{
    return std::get<0>(bar1s(ex, ep, ed, eq, 2));
}

ColVec bar1s(const RowVec &ex, const RowVec &ep, const RowVec &ed)
{
    return bar1s(ex, ep, ed, 0.0);
}

// ---- bar1we / bar1ws ----

std::tuple<Matrix, ColVec> bar1we(const RowVec &ex, const RowVec &ep, double eq)
{
    double DEA = ep(0) * ep(1);
    double kX = ep(2);
    double L = length1D(ex);

    Matrix K2(2, 2);
    K2 << 2.0, 1.0, //
        1.0, 2.0;

    Matrix Ke = axialStiffness(DEA, L) + kX * L / 6.0 * K2;
    return {Ke, axialLoad(eq, L)};
}

Matrix bar1we(const RowVec &ex, const RowVec &ep)
{
    return std::get<0>(bar1we(ex, ep, 0.0));
}

std::tuple<ColVec, ColVec, ColVec> bar1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep)
{
    double DEA = ep(0) * ep(1);
    auto s = axialSection(length1D(ex), ed(0), ed(1), DEA, eq, ep(2), nep, "bar1ws");
    return {s.N, s.u, s.X};
}

ColVec bar1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq)
{
    return std::get<0>(bar1ws(ex, ep, ed, eq, 2));
}

ColVec bar1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed)
{
    return bar1ws(ex, ep, ed, 0.0);
}

// ---- bar2e / bar2s ----

std::tuple<Matrix, ColVec> bar2e(const RowVec &ex, const RowVec &ey, const RowVec &ep, double eq)
{
    double DEA = ep(0) * ep(1);
    auto axis = axis2D(ex, ey);
    Matrix G = barTransform2D(axis);

    Matrix Ke = G.transpose() * axialStiffness(DEA, axis.L) * G;
    ColVec fe = G.transpose() * axialLoad(eq, axis.L);
    return {Ke, fe};
}

Matrix bar2e(const RowVec &ex, const RowVec &ey, const RowVec &ep)
{
    return std::get<0>(bar2e(ex, ey, ep, 0.0));
}

std::tuple<ColVec, ColVec, ColVec> bar2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                         double eq, int nep)
{
    double DEA = ep(0) * ep(1);
    auto axis = axis2D(ex, ey);

    ColVec a1 = barTransform2D(axis) * ed.transpose();
    auto s = axialSection(axis.L, a1(0), a1(1), DEA, eq, 0.0, nep, "bar2s");
    return {s.N, s.u, s.X};
}

ColVec bar2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, double eq)
{
    return std::get<0>(bar2s(ex, ey, ep, ed, eq, 2));
}

ColVec bar2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed)
{
    return bar2s(ex, ey, ep, ed, 0.0);
}

// ---- bar2ge / bar2gs ----

Matrix bar2ge(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX)
{
    double DEA = ep(0) * ep(1);
    auto axis = axis2D(ex, ey);
    double L = axis.L;

    Matrix K0le(4, 4);
    K0le << 1.0, 0.0, -1.0, 0.0, //
        0.0, 0.0, 0.0, 0.0,      //
        -1.0, 0.0, 1.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0;

    Matrix Ksle(4, 4);
    Ksle << 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, -1.0,    //
        0.0, 0.0, 0.0, 0.0,     //
        0.0, -1.0, 0.0, 1.0;

    Matrix Kle = DEA / L * K0le + QX / L * Ksle;
    Matrix G = nodeTransform2D(axis);

    return G.transpose() * Kle * G;
}

std::tuple<ColVec, double, ColVec, ColVec> bar2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep,
                                                  const RowVec &ed, int nep)
{
    double DEA = ep(0) * ep(1);
    auto axis = axis2D(ex, ey);

    ColVec edl = nodeTransform2D(axis) * ed.transpose();
    auto s = axialSection(axis.L, edl(0), edl(2), DEA, 0.0, 0.0, nep, "bar2gs");

    return {s.N, s.N(0), s.u, s.X};
}

std::tuple<ColVec, double> bar2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed)
{
    auto result = bar2gs(ex, ey, ep, ed, 2);
    return {std::get<0>(result), std::get<1>(result)};
}

} // namespace calfem
