//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Solid elements, ported from calfem-python (core.py).

#include <calfem/elements/solid.h>

#include "detail/gauss.h"
#include "detail/log.h"
#include "detail/shape_functions.h"

#include <cmath>
#include <limits>
#include <string>

namespace calfem {

using namespace detail;

namespace {

Matrix brickCoords(const RowVec &ex, const RowVec &ey, const RowVec &ez)
{
    Matrix c(8, 3);
    c.col(0) = ex.transpose();
    c.col(1) = ey.transpose();
    c.col(2) = ez.transpose();
    return c;
}

void checkJacobian(double detJ, const char *caller)
{
    if (detJ < 10.0 * std::numeric_limits<double>::epsilon())
        logMessage(std::string(caller) + ": warning, Jacobian determinant is zero or negative.");
}

// Strain-displacement matrix (6 x 24, strains [xx, yy, zz, xy, xz, yz]) and displacement
// interpolation (3 x 24)
void brickMatrices(const Matrix &dNx, const RowVec &N, Matrix &B, Matrix &N2)
{
    B = Matrix::Zero(6, 24);
    N2 = Matrix::Zero(3, 24);

    for (int k = 0; k < 8; k++)
    {
        B(0, 3 * k) = dNx(0, k);
        B(1, 3 * k + 1) = dNx(1, k);
        B(2, 3 * k + 2) = dNx(2, k);
        B(3, 3 * k) = dNx(1, k);
        B(3, 3 * k + 1) = dNx(0, k);
        B(4, 3 * k) = dNx(2, k);
        B(4, 3 * k + 2) = dNx(0, k);
        B(5, 3 * k + 1) = dNx(2, k);
        B(5, 3 * k + 2) = dNx(1, k);

        N2(0, 3 * k) = N(k);
        N2(1, 3 * k + 1) = N(k);
        N2(2, 3 * k + 2) = N(k);
    }
}

} // namespace

std::tuple<Matrix, ColVec> soli8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                  const Matrix &D, const RowVec &eq)
{
    auto rule = gaussRule3D(static_cast<int>(std::lround(ep(0))), "soli8e");
    Matrix coords = brickCoords(ex, ey, ez);
    ColVec b = eq.transpose();

    Matrix Ke = Matrix::Zero(24, 24);
    ColVec fe = ColVec::Zero(24);

    RowVec N;
    Matrix dNr, B, N2;

    for (Eigen::Index i = 0; i < rule.points.rows(); i++)
    {
        brick8Shape(rule.points(i, 0), rule.points(i, 1), rule.points(i, 2), N, dNr);

        Matrix JT = dNr * coords;
        double detJ = JT.determinant();
        checkJacobian(detJ, "soli8e");

        brickMatrices(JT.inverse() * dNr, N, B, N2);

        Ke += B.transpose() * D * B * detJ * rule.weights(i);
        fe += N2.transpose() * b * detJ * rule.weights(i);
    }

    return {Ke, fe};
}

Matrix soli8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(soli8e(ex, ey, ez, ep, D, RowVec::Zero(3)));
}

std::tuple<Matrix, Matrix, Matrix> soli8s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                          const Matrix &D, const RowVec &ed)
{
    auto rule = gaussRule3D(static_cast<int>(std::lround(ep(0))), "soli8s");
    Matrix coords = brickCoords(ex, ey, ez);
    auto ngp = rule.points.rows();

    Matrix et(ngp, 6), es(ngp, 6), eci(ngp, 3);

    RowVec N;
    Matrix dNr, B, N2;

    for (Eigen::Index i = 0; i < ngp; i++)
    {
        brick8Shape(rule.points(i, 0), rule.points(i, 1), rule.points(i, 2), N, dNr);
        eci.row(i) = N * coords;

        Matrix JT = dNr * coords;
        checkJacobian(JT.determinant(), "soli8s");

        brickMatrices(JT.inverse() * dNr, N, B, N2);

        ColVec ee = B * ed.transpose();
        et.row(i) = ee.transpose();
        es.row(i) = (D * ee).transpose();
    }

    return {et, es, eci};
}

} // namespace calfem
