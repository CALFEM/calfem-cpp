//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/elements/beam.h>

#include <cmath>

namespace calfem {

void beam3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep, const RowVec &eq,
            Matrix &Ke, ColVec &fe)
{
    ColVec bb(3);
    bb << ex(1) - ex(0), ey(1) - ey(0), ez(1) - ez(0);

    double L = sqrt((bb.transpose() * bb));

    RowVec n1 = (bb / L).transpose();

    double lc = sqrt((eo * eo.transpose()));
    RowVec n3 = eo / lc;

    double qx = eq(0);
    double qy = eq(1);
    double qz = eq(2);
    double qw = eq(3);

    double E = ep(0);
    double Gs = ep(1);
    double A = ep(2);
    double Iy = ep(3);
    double Iz = ep(4);
    double Kv = ep(5);

    double a = E * A / L;
    double b = 12.0 * E * Iz / pow(L, 3);
    double c = 6.0 * E * Iz / pow(L, 2);
    double d = 12.0 * E * Iy / pow(L, 3);
    double e = 6.0 * E * Iy / pow(L, 2);
    double f = Gs * Kv / L;
    double g = 2.0 * E * Iy / L;
    double h = 2 * E * Iz / L;

    Matrix Kle(12, 12);

    Kle << a, 0, 0, 0, 0, 0, -a, 0, 0, 0, 0, 0, 0, b, 0, 0, 0, c, 0, -b, 0, 0, 0, c, 0, 0, d, 0, -e, 0, 0, 0, -d, 0, -e,
        0, 0, 0, 0, f, 0, 0, 0, 0, 0, -f, 0, 0, 0, 0, -e, 0, 2 * g, 0, 0, 0, e, 0, g, 0, 0, c, 0, 0, 0, 2 * h, 0, -c, 0,
        0, 0, h, -a, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0, 0, -b, 0, 0, 0, -c, 0, b, 0, 0, 0, -c, 0, 0, -d, 0, e, 0, 0, 0, d,
        0, e, 0, 0, 0, 0, -f, 0, 0, 0, 0, 0, f, 0, 0, 0, 0, -e, 0, g, 0, 0, 0, e, 0, 2 * g, 0, 0, c, 0, 0, 0, h, 0, -c,
        0, 0, 0, 2 * h;

    ColVec fle(12);

    fle << qx, qy, qz, qw, -qz * L / 6.0, qy * L / 6.0, qx, qy, qz, qw, qz * L / 6.0, -qy * L / 6.0;

    fle = fle * L / 2.0;

    RowVec n2(3);

    n2 << n3(1) * n1(2) - n3(2) * n1(1), -n1(2) * n3(0) + n1(0) * n3(2), n3(0) * n1(1) - n1(0) * n3(1);

    // Use fixed-size matrix for stack allocation (eliminates heap overhead)
    Matrix3x3 An;

    An.row(0) = n1;
    An.row(1) = n2;
    An.row(2) = n3;

    // Use fixed-size matrix for stack allocation
    Matrix12x12 G = Matrix12x12::Zero();

    G.block(0, 0, 3, 3) = An;
    G.block(3, 3, 3, 3) = An;
    G.block(6, 6, 3, 3) = An;
    G.block(9, 9, 3, 3) = An;

    Ke = G.transpose() * Kle * G;
    fe = G.transpose() * fle;
}

void beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep, const RowVec &ed,
            const RowVec &eq, const int n, Matrix &es, Matrix &edi, ColVec &eci)
{
    ColVec bb(3);
    bb << ex(1) - ex(0), ey(1) - ey(0), ez(1) - ez(0);

    double L = sqrt((bb.transpose() * bb));

    RowVec n1 = (bb / L).transpose();

    double lc = sqrt((eo * eo.transpose()));
    RowVec n3 = eo / lc;

    double qx = eq(0);
    double qy = eq(1);
    double qz = eq(2);
    double qw = eq(3);

    double E = ep(0);
    double Gs = ep(1);
    double A = ep(2);
    double Iy = ep(3);
    double Iz = ep(4);
    double Kv = ep(5);

    RowVec n2(3);

    n2 << n3(1) * n1(2) - n3(2) * n1(1), -n1(2) * n3(0) + n1(0) * n3(2), n3(0) * n1(1) - n1(0) * n3(1);

    // Use fixed-size matrix for stack allocation (eliminates heap overhead)
    Matrix3x3 An;

    An.row(0) = n1;
    An.row(1) = n2;
    An.row(2) = n3;

    // Use fixed-size matrix for stack allocation
    Matrix12x12 G = Matrix12x12::Zero();

    G.block(0, 0, 3, 3) = An;
    G.block(3, 3, 3, 3) = An;
    G.block(6, 6, 3, 3) = An;
    G.block(9, 9, 3, 3) = An;

    // Use fixed-size vectors for stack allocation
    Vector12 u;
    Vector12 diffSol;

    diffSol << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -qx * pow(L, 2) / 2.0 / E / A, qy * pow(L, 4) / 24.0 / E / Iz,
        qz * pow(L, 4) / 24.0 / E / Iy, -qw * pow(L, 2) / 2.0 / Gs / Kv, -qz * pow(L, 3) / 6.0 / E / Iy,
        qy * pow(L, 3) / 6.0 / E / Iz;

    // u = G*ed.AsColumn() - diffSol;
    u = G * ed.transpose() - diffSol;

    // Use fixed-size matrix for stack allocation
    Matrix12x12 C;

    C << 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, L, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pow(L, 3),
        pow(L, 2), L, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pow(L, 3), pow(L, 2), L, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, L, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -3.0 * pow(L, 2), -2.0 * L, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0 * pow(L, 2), 2.0 * L, 1, 0.0, 0.0, 0.0, 0.0, 0.0,
        0, 0.0;

    Vector12 m;

    m = C.inverse() * u;

    es.resize(n, 6);
    es.setZero();

    edi.resize(n, 4);
    edi.setZero();

    eci.resize(n);
    eci.setZero();

    for (auto i = 0; i < n; i++)
    {
        eci(i) = double(i) * L / (n > 1 ? double(n - 1) : 1.0);
        double x = eci(i);

        // Use fixed-size matrices for stack allocation (critical for performance)
        Matrix6x12 T1;

        T1 << E * A, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0 * E * Iz, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0 * E * Iy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Gs * Kv, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0 * E * Iy * x,
            -2.0 * E * Iy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0 * E * Iz * x, 2.0 * E * Iz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0;

        Vector6 T2;

        T2 << -qx * x, -qy * x, -qz * x, -qw * x, -qz * pow(x, 2) / 2.0, qy * pow(x, 2) / 2.0;

        es.row(i) = (T1 * m + T2).transpose();

        Matrix4x12 T3;

        T3 << x, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pow(x, 3), pow(x, 2), x, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pow(x, 3), pow(x, 2), x, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, x, 1.0;

        Vector4 T4;

        T4 << -qx * pow(x, 2) / 2.0 / E / A, qy * pow(x, 4) / 24.0 / E / Iz, qz * pow(x, 4) / 24.0 / E / Iy,
            -qw * pow(x, 2) / 2.0 / Gs / Kv;

        edi.row(i) = (T3 * m + T4).transpose();
    }
}

// ---- calfem-python style ----
//
// Thin wrappers around the output-argument versions above, so both APIs give identical results.

Matrix beam3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep)
{
    Matrix Ke;
    ColVec fe;
    beam3e(ex, ey, ez, eo, ep, RowVec::Zero(4), Ke, fe);
    return Ke;
}

std::tuple<Matrix, ColVec> beam3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo,
                                  const RowVec &ep, const RowVec &eq)
{
    Matrix Ke;
    ColVec fe;
    beam3e(ex, ey, ez, eo, ep, eq, Ke, fe);
    return {Ke, fe};
}

Matrix beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep,
              const RowVec &ed)
{
    return beam3s(ex, ey, ez, eo, ep, ed, RowVec::Zero(4));
}

Matrix beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep,
              const RowVec &ed, const RowVec &eq)
{
    return std::get<0>(beam3s(ex, ey, ez, eo, ep, ed, eq, 2));
}

std::tuple<Matrix, Matrix, ColVec> beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo,
                                          const RowVec &ep, const RowVec &ed, const RowVec &eq, int nep)
{
    Matrix es, edi;
    ColVec eci;
    beam3s(ex, ey, ez, eo, ep, ed, eq, nep, es, edi, eci);
    return {es, edi, eci};
}

} // namespace calfem
