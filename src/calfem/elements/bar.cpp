//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/elements/bar.h>

#include <cmath>

namespace calfem {

void bar3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, double eq, Matrix &Ke, ColVec &fe)
{
    double E = ep(0);
    double A = ep(1);

    ColVec b(3);
    RowVec n(3);
    Matrix G(2, 6);
    Matrix Kle(2, 2);
    ColVec fle(2);

    b << ex(1) - ex(0), ey(1) - ey(0), ez(1) - ez(0);

    double L = sqrt((b.transpose() * b));

    n = b.transpose() / L;

    G << n(0), n(1), n(2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, n(0), n(1), n(2);

    double a = E * A / L;

    Kle << a, -a, -a, a;

    a = eq * L / 2.0;

    fle << a, a;

    Ke = G.transpose() * Kle * G;
    fe = G.transpose() * fle;
}

void bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const RowVec &ed, double eq, int n,
           ColVec &es, ColVec &edi, ColVec &eci)
{
    double E = ep(0);
    double A = ep(1);

    ColVec b(3);
    RowVec nn(3);
    Matrix G(2, 6);
    Matrix Kle(2, 2);
    ColVec fle(2);

    b << ex(1) - ex(0), ey(1) - ey(0), ez(1) - ez(0);

    double L = sqrt((b.transpose() * b));

    nn = b.transpose() / L;

    G << nn(0), nn(1), nn(2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nn(0), nn(1), nn(2);

    double c1 = E * A / L;

    Kle << c1, -c1, -c1, c1;

    c1 = eq * L / 2.0;

    fle << c1, c1;

    ColVec ae = ed.transpose();
    ColVec ale = G * ae;

    RowVec B(2);

    B << -1.0 / L, 1.0 / L;

    es.resize(n);
    edi.resize(n);
    eci.resize(n);

    double nIntervals = n > 1 ? n - 1.0 : 1.0;

    for (int i = 1; i <= n; i++)
    {
        eci(i - 1) = (i - 1.0) * L / nIntervals;
        double x = eci(i - 1);
        double up = -eq * (0.5 * pow(x, 2) - 0.5 * L * x) / E / A;
        double Np = -eq * (x - 0.5 * L);

        RowVec N(2);
        N << 1.0 - x / L, x / L;

        edi(i - 1) = N * ale + up; // [1 x 2] * [2 x 1]
        es(i - 1) = E * A * B * ale + Np;
        eci(i - 1) = x;
    }
}

// ---- calfem-python style ----
//
// Thin wrappers around the output-argument versions above, so both APIs give identical results.

Matrix bar3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep)
{
    Matrix Ke;
    ColVec fe;
    bar3e(ex, ey, ez, ep, 0.0, Ke, fe);
    return Ke;
}

std::tuple<Matrix, ColVec> bar3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, double eq)
{
    Matrix Ke;
    ColVec fe;
    bar3e(ex, ey, ez, ep, eq, Ke, fe);
    return {Ke, fe};
}

ColVec bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const RowVec &ed)
{
    return bar3s(ex, ey, ez, ep, ed, 0.0);
}

ColVec bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const RowVec &ed, double eq)
{
    return std::get<0>(bar3s(ex, ey, ez, ep, ed, eq, 2));
}

std::tuple<ColVec, ColVec, ColVec> bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                         const RowVec &ed, double eq, int nep)
{
    ColVec es, edi, eci;
    bar3s(ex, ey, ez, ep, ed, eq, nep, es, edi, eci);
    return {es, edi, eci};
}

} // namespace calfem
