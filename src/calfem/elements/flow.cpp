//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Field problem elements, ported from calfem-python (core.py).

#include <calfem/elements/flow.h>
#include <calfem/solvers.h>

#include "detail/gauss.h"
#include "detail/shape_functions.h"
#include "detail/log.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace calfem {

using namespace detail;

namespace {

// ---- Triangle ----

struct Triangle {
    Matrix B; // gradient matrix (2 x 3)
    double A; // area
};

Triangle triangle(const RowVec &ex, const RowVec &ey)
{
    Matrix C(3, 3);
    C << 1.0, ex(0), ey(0), //
        1.0, ex(1), ey(1),  //
        1.0, ex(2), ey(2);

    Matrix G(2, 3);
    G << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0;

    return {G * C.inverse(), 0.5 * C.determinant()};
}

// Evaluates the shape functions at Gauss point i of the rule
template <typename Shape> void evalShape(const Shape &shape, const GaussRule &rule, Eigen::Index i, RowVec &N, Matrix &dNr)
{
    if (rule.points.cols() == 2)
        shape(rule.points(i, 0), rule.points(i, 1), 0.0, N, dNr);
    else
        shape(rule.points(i, 0), rule.points(i, 1), rule.points(i, 2), N, dNr);
}

void checkJacobian(double detJ, const char *caller)
{
    if (detJ < 10.0 * std::numeric_limits<double>::epsilon())
        logMessage(std::string(caller) + ": warning, Jacobian determinant is zero or negative.");
}

// Stiffness matrix K = int B^T D B and load vector f = int N^T over the element.
template <typename Shape>
std::tuple<Matrix, ColVec> isoMatrices(const Matrix &coords, const Matrix &D, const GaussRule &rule,
                                       const Shape &shape, const char *caller)
{
    auto nnodes = coords.rows();
    Matrix Ke = Matrix::Zero(nnodes, nnodes);
    ColVec fe = ColVec::Zero(nnodes);

    RowVec N;
    Matrix dNr;

    for (Eigen::Index i = 0; i < rule.points.rows(); i++)
    {
        evalShape(shape, rule, i, N, dNr);

        Matrix JT = dNr * coords;
        double detJ = JT.determinant();
        checkJacobian(detJ, caller);

        Matrix B = JT.inverse() * dNr;
        Ke += B.transpose() * D * B * detJ * rule.weights(i);
        fe += N.transpose() * detJ * rule.weights(i);
    }

    return {Ke, fe};
}

// Flux es = -D B a, gradient et = B a and coordinates eci at the Gauss points.
template <typename Shape>
std::tuple<Matrix, Matrix, Matrix> isoFlux(const Matrix &coords, const Matrix &D, const GaussRule &rule,
                                           const Shape &shape, const Matrix &ed, const char *caller)
{
    auto ngp = rule.points.rows();
    auto dim = coords.cols();
    auto ncases = ed.rows();

    Matrix es(ngp * ncases, dim), et(ngp * ncases, dim), eci(ngp, dim);

    RowVec N;
    Matrix dNr;

    for (Eigen::Index i = 0; i < ngp; i++)
    {
        evalShape(shape, rule, i, N, dNr);
        eci.row(i) = N * coords;

        Matrix JT = dNr * coords;
        checkJacobian(JT.determinant(), caller);

        Matrix B = JT.inverse() * dNr;

        for (Eigen::Index j = 0; j < ncases; j++)
        {
            ColVec grad = B * ed.row(j).transpose();
            es.row(i + j * ngp) = (-D * grad).transpose();
            et.row(i + j * ngp) = grad.transpose();
        }
    }

    return {es, et, eci};
}

Matrix coords2D(const RowVec &ex, const RowVec &ey)
{
    Matrix c(ex.size(), 2);
    c.col(0) = ex.transpose();
    c.col(1) = ey.transpose();
    return c;
}

Matrix coords3D(const RowVec &ex, const RowVec &ey, const RowVec &ez)
{
    Matrix c(ex.size(), 3);
    c.col(0) = ex.transpose();
    c.col(1) = ey.transpose();
    c.col(2) = ez.transpose();
    return c;
}

int integrationRule(double ir)
{
    return static_cast<int>(std::lround(ir));
}

auto quad4 = [](double xsi, double eta, double, RowVec &N, Matrix &dNr) { quad4Shape(xsi, eta, N, dNr); };
auto quad8 = [](double xsi, double eta, double, RowVec &N, Matrix &dNr) { quad8Shape(xsi, eta, N, dNr); };
auto brick8 = [](double xsi, double eta, double zet, RowVec &N, Matrix &dNr) { brick8Shape(xsi, eta, zet, N, dNr); };

// ---- Quadrilateral from four triangles ----

// Triangles (local node numbers, node 4 is the centre)
const int quadTriangles[4][3] = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};

struct QuadTriangles {
    RowVec ex[4], ey[4];
};

QuadTriangles quadTriangleCoords(const RowVec &ex, const RowVec &ey)
{
    double xc = ex.sum() / 4.0, yc = ey.sum() / 4.0;
    RowVec x(5), y(5);
    x << ex(0), ex(1), ex(2), ex(3), xc;
    y << ey(0), ey(1), ey(2), ey(3), yc;

    QuadTriangles q;
    for (int t = 0; t < 4; t++)
    {
        q.ex[t].resize(3);
        q.ey[t].resize(3);
        for (int k = 0; k < 3; k++)
        {
            q.ex[t](k) = x(quadTriangles[t][k]);
            q.ey[t](k) = y(quadTriangles[t][k]);
        }
    }
    return q;
}

// Stiffness and load of the five-node (corner + centre) assembly of four triangles
std::tuple<Matrix, ColVec> quadAssembly(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                        double eq)
{
    auto q = quadTriangleCoords(ex, ey);

    Matrix K = Matrix::Zero(5, 5);
    ColVec f = ColVec::Zero(5);

    for (int t = 0; t < 4; t++)
    {
        auto [k1, f1] = flw2te(q.ex[t], q.ey[t], ep, D, eq);
        for (int i = 0; i < 3; i++)
        {
            f(quadTriangles[t][i]) += f1(i);
            for (int j = 0; j < 3; j++)
                K(quadTriangles[t][i], quadTriangles[t][j]) += k1(i, j);
        }
    }

    return {K, f};
}

} // namespace

// ---- flw2te / flw2ts ----

std::tuple<Matrix, ColVec> flw2te(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq)
{
    double t = ep(0);
    auto tri = triangle(ex, ey);

    Matrix Ke = tri.B.transpose() * D * tri.B * t * tri.A;
    ColVec fe = ColVec::Ones(3) * eq * tri.A * t / 3.0;

    return {Ke, fe};
}

Matrix flw2te(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(flw2te(ex, ey, ep, D, 0.0));
}

std::tuple<Matrix, Matrix> flw2ts(const Matrix &ex, const Matrix &ey, const Matrix &D, const Matrix &ed)
{
    Matrix es(ex.rows(), 2), et(ex.rows(), 2);

    for (Eigen::Index i = 0; i < ex.rows(); i++)
    {
        auto tri = triangle(ex.row(i), ey.row(i));
        ColVec grad = tri.B * ed.row(i).transpose();
        es.row(i) = (-D * grad).transpose();
        et.row(i) = grad.transpose();
    }

    return {es, et};
}

// ---- flw2qe / flw2qs ----

std::tuple<Matrix, ColVec> flw2qe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq)
{
    auto [K, f] = quadAssembly(ex, ey, ep, D, eq);

    IntColVec centre(1);
    centre << 5;
    return statcon(K, f, centre);
}

Matrix flw2qe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(flw2qe(ex, ey, ep, D, 0.0));
}

std::tuple<Matrix, Matrix> flw2qs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const Matrix &ed, double eq)
{
    auto [K, f] = quadAssembly(ex, ey, ep, D, eq);
    auto q = quadTriangleCoords(ex, ey);

    auto ncases = ed.rows();
    Matrix es = Matrix::Zero(ncases, 2), et = Matrix::Zero(ncases, 2);

    for (Eigen::Index j = 0; j < ncases; j++)
    {
        // Centre value from the condensed equation, corner values given
        ColVec a(5);
        a.head(4) = ed.row(j).transpose();
        a(4) = (f(4) - K.row(4).head(4).dot(a.head(4))) / K(4, 4);

        for (int t = 0; t < 4; t++)
        {
            RowVec at(3);
            for (int k = 0; k < 3; k++)
                at(k) = a(quadTriangles[t][k]);

            auto [s, g] = flw2ts(q.ex[t], q.ey[t], D, at);
            es.row(j) += s;
            et.row(j) += g;
        }
    }

    return {es / 4.0, et / 4.0};
}

std::tuple<Matrix, Matrix> flw2qs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const Matrix &ed)
{
    return flw2qs(ex, ey, ep, D, ed, 0.0);
}

// ---- flw2i4e / flw2i4s ----

std::tuple<Matrix, ColVec> flw2i4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq)
{
    double t = ep(0);
    auto rule = gaussRule2D(integrationRule(ep(1)), "flw2i4e");
    auto [Ke, fe] = isoMatrices(coords2D(ex, ey), D, rule, quad4, "flw2i4e");
    return {Ke * t, fe * t * eq};
}

Matrix flw2i4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(flw2i4e(ex, ey, ep, D, 0.0));
}

std::tuple<Matrix, Matrix, Matrix> flw2i4s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                           const Matrix &ed)
{
    auto rule = gaussRule2D(integrationRule(ep(1)), "flw2i4s");
    return isoFlux(coords2D(ex, ey), D, rule, quad4, ed, "flw2i4s");
}

// ---- flw2i8e / flw2i8s ----

std::tuple<Matrix, ColVec> flw2i8e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq)
{
    double t = ep(0);
    auto rule = gaussRule2D(integrationRule(ep(1)), "flw2i8e");
    auto [Ke, fe] = isoMatrices(coords2D(ex, ey), D, rule, quad8, "flw2i8e");
    return {Ke * t, fe * t * eq};
}

Matrix flw2i8e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(flw2i8e(ex, ey, ep, D, 0.0));
}

std::tuple<Matrix, Matrix, Matrix> flw2i8s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                           const Matrix &ed)
{
    auto rule = gaussRule2D(integrationRule(ep(1)), "flw2i8s");
    return isoFlux(coords2D(ex, ey), D, rule, quad8, ed, "flw2i8s");
}

// ---- flw3i8e / flw3i8s ----

std::tuple<Matrix, ColVec> flw3i8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                   const Matrix &D, double eq)
{
    auto rule = gaussRule3D(integrationRule(ep(0)), "flw3i8e");
    auto [Ke, fe] = isoMatrices(coords3D(ex, ey, ez), D, rule, brick8, "flw3i8e");
    return {Ke, fe * eq};
}

Matrix flw3i8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const Matrix &D)
{
    return std::get<0>(flw3i8e(ex, ey, ez, ep, D, 0.0));
}

std::tuple<Matrix, Matrix, Matrix> flw3i8s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                           const Matrix &D, const Matrix &ed)
{
    auto rule = gaussRule3D(integrationRule(ep(0)), "flw3i8s");
    return isoFlux(coords3D(ex, ey, ez), D, rule, brick8, ed, "flw3i8s");
}

} // namespace calfem
