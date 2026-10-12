//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Plate elements, ported from calfem-python (core.py).

#include <calfem/elements/plate.h>

#include <initializer_list>

namespace calfem {

std::tuple<Matrix, ColVec> platre(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq)
{
    double Lx = ex(2) - ex(0);
    double Ly = ey(2) - ey(0);
    double t = ep(0);

    Matrix Dp = t * t * t / 12.0 * D;
    double D11 = Dp(0, 0), D12 = Dp(0, 1), D22 = Dp(1, 1), D33 = Dp(2, 2);

    double A1 = Ly / (Lx * Lx * Lx);
    double A2 = Lx / (Ly * Ly * Ly);
    double A3 = 1.0 / Lx / Ly;
    double A4 = Ly / (Lx * Lx);
    double A5 = Lx / (Ly * Ly);
    double A6 = 1.0 / Lx;
    double A7 = 1.0 / Ly;
    double A8 = Ly / Lx;
    double A9 = Lx / Ly;

    double C1 = 4 * A1 * D11 + 4 * A2 * D22 + 2 * A3 * D12 + 5.6 * A3 * D33;
    double C2 = -4 * A1 * D11 + 2 * A2 * D22 - 2 * A3 * D12 - 5.6 * A3 * D33;
    double C3 = 2 * A1 * D11 - 4 * A2 * D22 - 2 * A3 * D12 - 5.6 * A3 * D33;
    double C4 = -2 * A1 * D11 - 2 * A2 * D22 + 2 * A3 * D12 + 5.6 * A3 * D33;
    double C5 = 2 * A5 * D22 + A6 * D12 + 0.4 * A6 * D33;
    double C6 = 2 * A4 * D11 + A7 * D12 + 0.4 * A7 * D33;
    double C7 = 2 * A5 * D22 + 0.4 * A6 * D33;
    double C8 = 2 * A4 * D11 + 0.4 * A7 * D33;
    double C9 = A5 * D22 - A6 * D12 - 0.4 * A6 * D33;
    double C10 = A4 * D11 - A7 * D12 - 0.4 * A7 * D33;
    double C11 = A5 * D22 - 0.4 * A6 * D33;
    double C12 = A4 * D11 - 0.4 * A7 * D33;
    double C13 = 4.0 / 3.0 * A9 * D22 + 8.0 / 15.0 * A8 * D33;
    double C14 = 4.0 / 3.0 * A8 * D11 + 8.0 / 15.0 * A9 * D33;
    double C15 = 2.0 / 3.0 * A9 * D22 - 8.0 / 15.0 * A8 * D33;
    double C16 = 2.0 / 3.0 * A8 * D11 - 8.0 / 15.0 * A9 * D33;
    double C17 = 2.0 / 3.0 * A9 * D22 - 2.0 / 15.0 * A8 * D33;
    double C18 = 2.0 / 3.0 * A8 * D11 - 2.0 / 15.0 * A9 * D33;
    double C19 = 1.0 / 3.0 * A9 * D22 + 2.0 / 15.0 * A8 * D33;
    double C20 = 1.0 / 3.0 * A8 * D11 + 2.0 / 15.0 * A9 * D33;
    double C21 = D12;

    // Upper triangle, row by row from the diagonal
    const std::initializer_list<double> upper[12] = {
        {C1, C5, -C6, C2, C9, -C8, C4, C11, -C12, C3, C7, -C10},
        {C13, -C21, C9, C15, 0, -C11, C19, 0, -C7, C17, 0},
        {C14, C8, 0, C18, C12, 0, C20, -C10, 0, C16},
        {C1, C5, C6, C3, C7, C10, C4, C11, C12},
        {C13, C21, -C7, C17, 0, -C11, C19, 0},
        {C14, C10, 0, C16, -C12, 0, C20},
        {C1, -C5, C6, C2, -C9, C8},
        {C13, -C21, -C9, C15, 0},
        {C14, -C8, 0, C18},
        {C1, -C5, -C6},
        {C13, C21},
        {C14}};

    Matrix Ke = Matrix::Zero(12, 12);
    for (int i = 0; i < 12; i++)
    {
        int j = i;
        for (double value : upper[i])
        {
            Ke(i, j) = value;
            Ke(j, i) = value;
            j++;
        }
    }

    double R1 = eq * Lx * Ly / 4.0;
    double R2 = eq * Lx * Ly * Ly / 24.0;
    double R3 = eq * Ly * Lx * Lx / 24.0;

    ColVec fe(12);
    fe << R1, R2, -R3, R1, R2, R3, R1, -R2, R3, R1, -R2, -R3;

    return {Ke, fe};
}

Matrix platre(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(platre(ex, ey, ep, D, 0.0));
}

} // namespace calfem
