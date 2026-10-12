//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Curve evaluation and discretization for the built-in mesher (not installed).

#pragma once

#include <calfem/geometry.h>

#include <Eigen/Dense>

#include <vector>

namespace calfem::detail {

using Vec3 = Eigen::Vector3d;

// Parametric shape of a geometry curve, x(t) for t in [0, 1], with arc length lookup.
class CurveShape {
public:
    CurveShape(const Geometry &geometry, int curveID);

    Vec3 eval(double t) const;
    double length() const;

    // Point at the arc length fraction s in [0, 1]
    Vec3 atLength(double s) const;

private:
    Vec3 evalSpline(double t) const;
    Vec3 evalBSpline(double t) const;

    Geometry::CurveType m_type;
    std::vector<Vec3> m_points; // defining points (spline/bspline control points)
    bool m_closed = false;

    // Circle/ellipse arc: centre, axes and angles
    Vec3 m_centre, m_e1, m_e2;
    double m_a = 0.0, m_b = 0.0, m_t0 = 0.0, m_t1 = 0.0;

    // Dense samples for arc length lookup
    std::vector<double> m_t, m_s;
};

// Arc length fractions xi_0 = 0 < xi_1 < ... < xi_n = 1 of the nodes on a curve.
//
// With el_on_curve > 0 the curve gets that many elements, distributed uniformly or according to
// the curve's Progression/Bump setting. Otherwise the number of elements follows from the
// element sizes h0, h1 at the curve ends (linearly varying along the curve). With even, the
// number of elements computed from the sizes is rounded up to an even number.
std::vector<double> curveDivision(const Geometry::Curve &curve, double length, double h0, double h1,
                                  bool even = false);

} // namespace calfem::detail
