//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Plane stress/strain elements, ported from calfem-python (core.py).

#include <calfem/elements/plane.h>
#include <calfem/solvers.h>

#include "detail/gauss.h"
#include "detail/log.h"
#include "detail/shape_functions.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace calfem {

using namespace detail;

namespace {

// In-plane components [xx, yy, xy] of a 4-component vector [xx, yy, zz, xy]
const std::vector<Eigen::Index> planeComponents = {0, 1, 3};

int planeType(const RowVec &ep, const char *caller)
{
    auto ptype = static_cast<int>(std::lround(ep(0)));
    if ((ptype != 1) && (ptype != 2))
        throw std::invalid_argument(std::string(caller) + ": ptype must be 1 (plane stress) or 2 (plane strain).");
    return ptype;
}

// 3 x 3 in-plane material matrix. A 4 x 4 D is condensed for plane stress (zz stress zero)
// and reduced for plane strain (zz strain zero).
Matrix planeMaterial(int ptype, const Matrix &D)
{
    if (D.cols() <= 3)
        return D;

    if (ptype == 1)
    {
        Matrix Cm = D.inverse();
        return Matrix(Cm(planeComponents, planeComponents)).inverse();
    }

    return D(planeComponents, planeComponents);
}

// Stresses and strains with as many components as D from the in-plane strains e = [xx, yy, xy]
void planeStressStrain(int ptype, const Matrix &D, const ColVec &e, ColVec &ss, ColVec &ee)
{
    if (D.cols() <= 3)
    {
        ee = e;
        ss = D * e;
    }
    else if (ptype == 1)
    {
        // Plane stress: zz stress zero, zz strain from the compliance
        Matrix Cm = D.inverse();
        Matrix Dm = Matrix(Cm(planeComponents, planeComponents)).inverse();
        ss = ColVec::Zero(D.cols());
        ss(planeComponents) = Dm * e;
        ee = Cm * ss;
    }
    else
    {
        // Plane strain: zz strain zero
        ee = ColVec::Zero(D.cols());
        ee(planeComponents) = e;
        ss = D * ee;
    }
}

// ---- Constant strain triangle ----

struct PlaneTriangle {
    Matrix B; // 3 x 6
    double A; // area
};

PlaneTriangle planeTriangle(const RowVec &ex, const RowVec &ey)
{
    Matrix C(6, 6);
    C << 1, ex(0), ey(0), 0, 0, 0, //
        0, 0, 0, 1, ex(0), ey(0),  //
        1, ex(1), ey(1), 0, 0, 0,  //
        0, 0, 0, 1, ex(1), ey(1),  //
        1, ex(2), ey(2), 0, 0, 0,  //
        0, 0, 0, 1, ex(2), ey(2);

    Matrix S(3, 6);
    S << 0, 1, 0, 0, 0, 0, //
        0, 0, 0, 0, 0, 1,  //
        0, 0, 1, 0, 1, 0;

    Matrix Ca(3, 3);
    Ca << 1, ex(0), ey(0), //
        1, ex(1), ey(1),   //
        1, ex(2), ey(2);

    return {S * C.inverse(), 0.5 * Ca.determinant()};
}

const RowVec &zeroLoad2()
{
    static const RowVec zero = RowVec::Zero(2);
    return zero;
}

// ---- Quadrilateral from four triangles (corner dofs 0..7, centre dofs 8, 9) ----

const int quadTriangleDofs[4][6] = {{0, 1, 2, 3, 8, 9}, {2, 3, 4, 5, 8, 9}, {4, 5, 6, 7, 8, 9}, {6, 7, 0, 1, 8, 9}};
const int quadTriangleNodes[4][3] = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};

struct QuadTriangles {
    RowVec ex[4], ey[4];
};

QuadTriangles quadTriangles(const RowVec &ex, const RowVec &ey)
{
    RowVec x(5), y(5);
    x << ex(0), ex(1), ex(2), ex(3), ex.sum() / 4.0;
    y << ey(0), ey(1), ey(2), ey(3), ey.sum() / 4.0;

    QuadTriangles q;
    for (int t = 0; t < 4; t++)
    {
        q.ex[t].resize(3);
        q.ey[t].resize(3);
        for (int k = 0; k < 3; k++)
        {
            q.ex[t](k) = x(quadTriangleNodes[t][k]);
            q.ey[t](k) = y(quadTriangleNodes[t][k]);
        }
    }
    return q;
}

std::tuple<Matrix, ColVec> quadAssembly(const QuadTriangles &q, const RowVec &ep, const Matrix &D, const RowVec &eq)
{
    Matrix K = Matrix::Zero(10, 10);
    ColVec f = ColVec::Zero(10);

    for (int t = 0; t < 4; t++)
    {
        auto [k1, f1] = plante(q.ex[t], q.ey[t], ep, D, eq);
        for (int i = 0; i < 6; i++)
        {
            f(quadTriangleDofs[t][i]) += f1(i);
            for (int j = 0; j < 6; j++)
                K(quadTriangleDofs[t][i], quadTriangleDofs[t][j]) += k1(i, j);
        }
    }

    return {K, f};
}

void checkJacobian(double detJ, const char *caller)
{
    if (detJ < 10.0 * std::numeric_limits<double>::epsilon())
        logMessage(std::string(caller) + ": warning, Jacobian determinant is zero or negative.");
}

// Strain-displacement matrix (3 x 8) and displacement interpolation (2 x 8) of the 4-node quad
void quad4Matrices(const Matrix &dNx, const RowVec &N, Matrix &B, Matrix &N2)
{
    B = Matrix::Zero(3, 8);
    N2 = Matrix::Zero(2, 8);

    for (int k = 0; k < 4; k++)
    {
        B(0, 2 * k) = dNx(0, k);
        B(1, 2 * k + 1) = dNx(1, k);
        B(2, 2 * k) = dNx(1, k);
        B(2, 2 * k + 1) = dNx(0, k);
        N2(0, 2 * k) = N(k);
        N2(1, 2 * k + 1) = N(k);
    }
}

Matrix quadCoords(const RowVec &ex, const RowVec &ey)
{
    Matrix c(4, 2);
    c.col(0) = ex.transpose();
    c.col(1) = ey.transpose();
    return c;
}

} // namespace

// ---- plante / plants / plantf ----

std::tuple<Matrix, ColVec> plante(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &eq)
{
    int ptype = planeType(ep, "plante");
    double t = ep(1);
    double bx = eq(0), by = eq(1);

    auto tri = planeTriangle(ex, ey);
    Matrix Dm = planeMaterial(ptype, D);

    Matrix Ke = tri.B.transpose() * Dm * tri.B * tri.A * t;

    ColVec fe(6);
    fe << bx, by, bx, by, bx, by;
    fe *= tri.A / 3.0 * t;

    return {Ke, fe};
}

Matrix plante(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(plante(ex, ey, ep, D, zeroLoad2()));
}

std::tuple<Matrix, Matrix> plants(const Matrix &ex, const Matrix &ey, const RowVec &ep, const Matrix &D,
                                  const Matrix &ed)
{
    int ptype = planeType(ep, "plants");

    if ((ex.rows() != 1) && (ex.rows() != ed.rows()))
        throw std::invalid_argument("plants: ex must have one row, or one row per row of ed.");

    Matrix es(ed.rows(), D.cols()), et(ed.rows(), D.cols());

    for (Eigen::Index i = 0; i < ed.rows(); i++)
    {
        auto ie = ex.rows() == 1 ? 0 : i;
        auto tri = planeTriangle(ex.row(ie), ey.row(ie));

        ColVec ss, ee;
        planeStressStrain(ptype, D, tri.B * ed.row(i).transpose(), ss, ee);

        es.row(i) = ss.transpose();
        et.row(i) = ee.transpose();
    }

    return {es, et};
}

ColVec plantf(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &es)
{
    planeType(ep, "plantf");
    double t = ep(1);

    ColVec stress(3);
    if (es.size() == 3)
        stress = es.transpose();
    else if (es.size() == 4)
        stress << es(0), es(1), es(3);
    else
        throw std::invalid_argument("plantf: es must contain 3 or 4 stress components.");

    auto tri = planeTriangle(ex, ey);
    return tri.A * t * (tri.B.transpose() * stress);
}

// ---- planqe / planqs ----

std::tuple<Matrix, ColVec> planqe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &eq)
{
    planeType(ep, "planqe");
    auto [K, f] = quadAssembly(quadTriangles(ex, ey), ep, D, eq);

    IntColVec centre(2);
    centre << 9, 10;
    return statcon(K, f, centre);
}

Matrix planqe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(planqe(ex, ey, ep, D, zeroLoad2()));
}

std::tuple<RowVec, RowVec> planqs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &ed, const RowVec &eq)
{
    planeType(ep, "planqs");

    if ((ex.size() != 4) || (ey.size() != 4) || (ed.size() != 8))
        throw std::invalid_argument("planqs: one element at a time (ex, ey with 4 and ed with 8 values).");

    auto q = quadTriangles(ex, ey);
    auto [K, f] = quadAssembly(q, ep, D, eq);

    // Centre displacements from the condensed equations, corner displacements given
    ColVec a(10);
    a.head(8) = ed.transpose();
    a.tail(2) = K.bottomRightCorner(2, 2).ldlt().solve(ColVec(f.tail(2) - K.bottomLeftCorner(2, 8) * a.head(8)));

    RowVec es = RowVec::Zero(D.cols()), et = RowVec::Zero(D.cols());
    double Atot = 0.0;

    for (int t = 0; t < 4; t++)
    {
        RowVec edt(6);
        for (int k = 0; k < 6; k++)
            edt(k) = a(quadTriangleDofs[t][k]);

        double A = planeTriangle(q.ex[t], q.ey[t]).A;
        auto [s, e] = plants(q.ex[t], q.ey[t], ep, D, edt);

        es += A * s.row(0);
        et += A * e.row(0);
        Atot += A;
    }

    return {es / Atot, et / Atot};
}

std::tuple<RowVec, RowVec> planqs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &ed)
{
    return planqs(ex, ey, ep, D, ed, zeroLoad2());
}

// ---- plani4e / plani4s ----

std::tuple<Matrix, ColVec> plani4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                   const RowVec &eq)
{
    int ptype = planeType(ep, "plani4e");
    double t = ep(1);
    auto rule = gaussRule2D(static_cast<int>(std::lround(ep(2))), "plani4e");

    Matrix Dm = planeMaterial(ptype, D);
    Matrix coords = quadCoords(ex, ey);
    ColVec q = eq.transpose();

    Matrix Ke = Matrix::Zero(8, 8);
    ColVec fe = ColVec::Zero(8);

    RowVec N;
    Matrix dNr, B, N2;

    for (Eigen::Index i = 0; i < rule.points.rows(); i++)
    {
        quad4Shape(rule.points(i, 0), rule.points(i, 1), N, dNr);

        Matrix JT = dNr * coords;
        double detJ = JT.determinant();
        checkJacobian(detJ, "plani4e");

        quad4Matrices(JT.inverse() * dNr, N, B, N2);

        Ke += B.transpose() * Dm * B * detJ * rule.weights(i) * t;
        fe += N2.transpose() * q * detJ * rule.weights(i) * t;
    }

    return {Ke, fe};
}

Matrix plani4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(plani4e(ex, ey, ep, D, zeroLoad2()));
}

std::tuple<Matrix, Matrix> plani4s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                   const RowVec &ed)
{
    int ptype = planeType(ep, "plani4s");
    auto rule = gaussRule2D(static_cast<int>(std::lround(ep(2))), "plani4s");

    Matrix coords = quadCoords(ex, ey);
    auto ngp = rule.points.rows();

    Matrix es(ngp, D.cols()), et(ngp, D.cols());

    RowVec N;
    Matrix dNr, B, N2;

    for (Eigen::Index i = 0; i < ngp; i++)
    {
        quad4Shape(rule.points(i, 0), rule.points(i, 1), N, dNr);

        Matrix JT = dNr * coords;
        checkJacobian(JT.determinant(), "plani4s");

        quad4Matrices(JT.partialPivLu().solve(dNr), N, B, N2);

        ColVec ss, ee;
        planeStressStrain(ptype, D, B * ed.transpose(), ss, ee);

        es.row(i) = ss.transpose();
        et.row(i) = ee.transpose();
    }

    return {es, et};
}

} // namespace calfem
