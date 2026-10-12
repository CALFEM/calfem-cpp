//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Internal helpers shared by the element routines (not installed).

#pragma once

#include <calfem/types.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace calfem::detail {

// Evaluation points x_i = i * (L / (ne - 1)), i = 0..ne-1, along an element of length L.
// Same values as calfem-python's np.arange(0, L + L/(ne-1), L/(ne-1)).

inline ColVec evalPoints(double L, int ne, const char *caller)
{
    if (ne < 2)
        throw std::invalid_argument(std::string(caller) + ": nep must be at least 2.");

    ColVec X(ne);
    double step = L / (ne - 1);

    for (int i = 0; i < ne; i++)
        X(i) = i * step;

    return X;
}

// Length and direction cosines of a 2D element.

struct Axis2D {
    double L;
    double nxX; // cos of the angle between the local x-axis and the global X-axis
    double nyX; // cos of the angle between the local x-axis and the global Y-axis
};

inline Axis2D axis2D(const RowVec &ex, const RowVec &ey)
{
    double dx = ex(1) - ex(0);
    double dy = ey(1) - ey(0);
    double L = std::sqrt(dx * dx + dy * dy);
    return {L, dx / L, dy / L};
}

// Global-to-local transformation for a 2D bar with axial dofs only (2 x 4).

inline Matrix barTransform2D(const Axis2D &a)
{
    Matrix G(2, 4);
    G << a.nxX, a.nyX, 0.0, 0.0, //
        0.0, 0.0, a.nxX, a.nyX;
    return G;
}

// Global-to-local transformation for two nodes with (ux, uy) dofs (4 x 4).

inline Matrix nodeTransform2D(const Axis2D &a)
{
    double nxY = -a.nyX;
    double nyY = a.nxX;

    Matrix G(4, 4);
    G << a.nxX, a.nyX, 0.0, 0.0, //
        nxY, nyY, 0.0, 0.0,      //
        0.0, 0.0, a.nxX, a.nyX,  //
        0.0, 0.0, nxY, nyY;
    return G;
}

// Global-to-local transformation for a 2D beam with (ux, uy, theta) dofs (6 x 6).

inline Matrix beamTransform2D(const Axis2D &a)
{
    double nxY = -a.nyX;
    double nyY = a.nxX;

    Matrix G(6, 6);
    G << a.nxX, a.nyX, 0.0, 0.0, 0.0, 0.0, //
        nxY, nyY, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 0.0, a.nxX, a.nyX, 0.0,  //
        0.0, 0.0, 0.0, nxY, nyY, 0.0,      //
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    return G;
}

// Coefficients of the linear axial displacement u(x) = c0 + c1 x from the end values.

inline Matrix axialInterpolation(double L)
{
    Matrix C1(2, 2);
    C1 << 1.0, 0.0, //
        -1.0 / L, 1.0 / L;
    return C1;
}

// Coefficients of the cubic deflection v(x) = c0 + c1 x + c2 x^2 + c3 x^3 from
// [v1, theta1, v2, theta2].

inline Matrix bendingInterpolation(double L)
{
    Matrix C2(4, 4);
    C2 << 1.0, 0.0, 0.0, 0.0,                                                //
        0.0, 1.0, 0.0, 0.0,                                                  //
        -3.0 / (L * L), -2.0 / L, 3.0 / (L * L), -1.0 / L,                   //
        2.0 / (L * L * L), 1.0 / (L * L), -2.0 / (L * L * L), 1.0 / (L * L);
    return C2;
}

struct AxialSection {
    ColVec N; // normal force
    ColVec u; // axial displacement
    ColVec X; // evaluation points
};

// Normal force and axial displacement along a bar with end displacements u1, u2, distributed
// load qX and (optionally) an elastic foundation kX.

inline AxialSection axialSection(double L, double u1, double u2, double DEA, double qX, double kX, int ne,
                          const char *caller)
{
    ColVec X = evalPoints(L, ne, caller);

    ColVec a1(2);
    a1 << u1, u2;
    ColVec C1a = axialInterpolation(L) * a1;

    Eigen::ArrayXd x = X.array();
    Eigen::ArrayXd u = C1a(0) + x * C1a(1);
    Eigen::ArrayXd du = Eigen::ArrayXd::Constant(ne, C1a(1));

    if (DEA != 0.0)
    {
        if (kX != 0.0)
        {
            u += kX / DEA * ((x.square() - L * x) / 2.0 * C1a(0) + (x.cube() - L * L * x) / 6.0 * C1a(1));
            du += kX / DEA * ((2.0 * x - L) / 2.0 * C1a(0) + (3.0 * x.square() - L * L) / 6.0 * C1a(1));
        }

        u -= (x.square() - L * x) * qX / (2.0 * DEA);
        du -= (2.0 * x - L) * qX / (2.0 * DEA);
    }

    return {ColVec(DEA * du), ColVec(u), X};
}

struct BendingSection {
    Eigen::ArrayXd v;   // deflection
    Eigen::ArrayXd dv;  // slope
    Eigen::ArrayXd d2v; // curvature
    Eigen::ArrayXd d3v; // third derivative
};

// Deflection and derivatives along a beam from the end values a2 = [v1, theta1, v2, theta2],
// with distributed load qY and (optionally) an elastic foundation kY, at the points X.

inline BendingSection bendingSection(double L, const ColVec &a2, double DEI, double qY, double kY, const ColVec &X)
{
    ColVec C2a = bendingInterpolation(L) * a2;
    double c0 = C2a(0), c1 = C2a(1), c2 = C2a(2), c3 = C2a(3);

    Eigen::ArrayXd x = X.array();
    Eigen::ArrayXd x2 = x.square();
    Eigen::ArrayXd x3 = x2 * x;
    double L2 = L * L, L3 = L2 * L, L4 = L3 * L, L5 = L4 * L;

    BendingSection b;
    b.v = c0 + c1 * x + c2 * x2 + c3 * x3;
    b.dv = c1 + 2.0 * c2 * x + 3.0 * c3 * x2;
    b.d2v = Eigen::ArrayXd::Constant(x.size(), 2.0 * c2) + 6.0 * c3 * x;
    b.d3v = Eigen::ArrayXd::Constant(x.size(), 6.0 * c3);

    if (DEI != 0.0)
    {
        if (kY != 0.0)
        {
            Eigen::ArrayXd x4 = x3 * x, x5 = x4 * x, x6 = x5 * x, x7 = x6 * x;

            b.v -= kY / DEI *
                   ((x4 - 2.0 * L * x3 + L2 * x2) / 24.0 * c0 + (x5 - 3.0 * L2 * x3 + 2.0 * L3 * x2) / 120.0 * c1 +
                    (x6 - 4.0 * L3 * x3 + 3.0 * L4 * x2) / 360.0 * c2 +
                    (x7 - 5.0 * L4 * x3 + 4.0 * L5 * x2) / 840.0 * c3);
            b.d2v -= kY / DEI *
                     ((6.0 * x2 - 6.0 * L * x + L2) / 12.0 * c0 + (10.0 * x3 - 9.0 * L2 * x + 2.0 * L3) / 60.0 * c1 +
                      (5.0 * x4 - 4.0 * L3 * x + L4) / 60.0 * c2 +
                      (21.0 * x5 - 15.0 * L4 * x + 4.0 * L5) / 420.0 * c3);
            b.d3v -= kY / DEI *
                     ((2.0 * x - L) / 2.0 * c0 + (10.0 * x2 - 3.0 * L2) / 20.0 * c1 + (5.0 * x3 - L3) / 15.0 * c2 +
                      (7.0 * x4 - L4) / 28.0 * c3);
        }

        b.v += (x2 * x2 - 2.0 * L * x3 + L2 * x2) * qY / (24.0 * DEI);
        b.dv += (2.0 * x3 - 3.0 * L * x2 + L2 * x) * qY / (12.0 * DEI);
        b.d2v += (6.0 * x2 - 6.0 * L * x + L2) * qY / (12.0 * DEI);
        b.d3v += (2.0 * x - L) * qY / (2.0 * DEI);
    }

    return b;
}

} // namespace calfem::detail
