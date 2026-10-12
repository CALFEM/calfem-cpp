//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// 1D beam elements, ported from calfem-python (core.py).

#include <calfem/elements/beam.h>

#include "detail/element_utils.h"

#include <cmath>

namespace calfem {

using namespace detail;

namespace {

Matrix bendingStiffness1D(double DEI, double L)
{
    Matrix K(4, 4);
    K << 12.0, 6.0 * L, -12.0, 6.0 * L,               //
        6.0 * L, 4.0 * L * L, -6.0 * L, 2.0 * L * L,  //
        -12.0, -6.0 * L, 12.0, -6.0 * L,              //
        6.0 * L, 2.0 * L * L, -6.0 * L, 4.0 * L * L;
    return DEI / (L * L * L) * K;
}

ColVec bendingLoad1D(double qY, double L)
{
    ColVec f(4);
    f << L / 2.0, L * L / 12.0, L / 2.0, -L * L / 12.0;
    return qY * f;
}

std::tuple<Matrix, ColVec, ColVec> beam1Section(const RowVec &ex, double DEI, const RowVec &ed, double qY, double kY,
                                                int nep, const char *caller)
{
    double L = std::abs(ex(1) - ex(0));
    ColVec X = evalPoints(L, nep, caller);

    auto b = bendingSection(L, ed.transpose(), DEI, qY, kY, X);

    Matrix es(nep, 2);
    es.col(0) = -DEI * b.d3v; // V
    es.col(1) = DEI * b.d2v;  // M

    return {es, ColVec(b.v), X};
}

} // namespace

// ---- beam1e / beam1s ----

std::tuple<Matrix, ColVec> beam1e(const RowVec &ex, const RowVec &ep, double eq)
{
    double DEI = ep(0) * ep(1);
    double L = std::abs(ex(1) - ex(0));
    return {bendingStiffness1D(DEI, L), bendingLoad1D(eq, L)};
}

Matrix beam1e(const RowVec &ex, const RowVec &ep)
{
    return std::get<0>(beam1e(ex, ep, 0.0));
}

std::tuple<Matrix, ColVec, ColVec> beam1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep)
{
    return beam1Section(ex, ep(0) * ep(1), ed, eq, 0.0, nep, "beam1s");
}

Matrix beam1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq)
{
    return std::get<0>(beam1s(ex, ep, ed, eq, 2));
}

Matrix beam1s(const RowVec &ex, const RowVec &ep, const RowVec &ed)
{
    return beam1s(ex, ep, ed, 0.0);
}

// ---- beam1we / beam1ws ----

std::tuple<Matrix, ColVec> beam1we(const RowVec &ex, const RowVec &ep, double eq)
{
    double DEI = ep(0) * ep(1);
    double kY = ep(2);
    double L = std::abs(ex(1) - ex(0));

    Matrix Ks(4, 4);
    Ks << 156.0, 22.0 * L, 54.0, -13.0 * L,                //
        22.0 * L, 4.0 * L * L, 13.0 * L, -3.0 * L * L,     //
        54.0, 13.0 * L, 156.0, -22.0 * L,                  //
        -13.0 * L, -3.0 * L * L, -22.0 * L, 4.0 * L * L;

    Matrix Ke = bendingStiffness1D(DEI, L) + kY * L / 420.0 * Ks;
    return {Ke, bendingLoad1D(eq, L)};
}

Matrix beam1we(const RowVec &ex, const RowVec &ep)
{
    return std::get<0>(beam1we(ex, ep, 0.0));
}

std::tuple<Matrix, ColVec, ColVec> beam1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep)
{
    return beam1Section(ex, ep(0) * ep(1), ed, eq, ep(2), nep, "beam1ws");
}

Matrix beam1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq)
{
    return std::get<0>(beam1ws(ex, ep, ed, eq, 2));
}

Matrix beam1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed)
{
    return beam1ws(ex, ep, ed, 0.0);
}

} // namespace calfem
