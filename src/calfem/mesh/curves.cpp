//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include "curves.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace calfem::detail {

namespace {

const double pi = 3.14159265358979323846;
const int samplesPerCurve = 512;

Vec3 pointCoord(const Geometry &geometry, int ID, int curveID)
{
    auto it = geometry.points().find(ID);
    if (it == geometry.points().end())
        throw std::invalid_argument("Mesh: curve " + std::to_string(curveID) + " uses point " + std::to_string(ID) +
                                    ", which does not exist.");
    return it->second.coord.transpose();
}

// Positions where the integral of w(xi) over [0, xi_k] is k/n of the total (w > 0)
template <typename W> std::vector<double> equalIntegral(int n, const W &w)
{
    const int m = 2000;
    std::vector<double> cum(m + 1, 0.0);
    for (int i = 1; i <= m; i++)
    {
        double a = double(i - 1) / m, b = double(i) / m;
        cum[i] = cum[i - 1] + 0.5 * (w(a) + w(b)) / m;
    }

    std::vector<double> xi(n + 1);
    xi[0] = 0.0;
    xi[n] = 1.0;
    for (int k = 1; k < n; k++)
    {
        double target = cum[m] * k / n;
        auto it = std::lower_bound(cum.begin(), cum.end(), target);
        auto i = static_cast<int>(it - cum.begin());
        double f = (target - cum[i - 1]) / (cum[i] - cum[i - 1]);
        xi[k] = (i - 1 + f) / m;
    }
    return xi;
}

} // namespace

CurveShape::CurveShape(const Geometry &geometry, int curveID)
{
    const auto &curve = geometry.curves().at(curveID);
    m_type = curve.type;

    for (int ID : curve.points)
        m_points.push_back(pointCoord(geometry, ID, curveID));

    m_closed = (curve.points.size() > 2) && (curve.points.front() == curve.points.back());

    if (m_type == Geometry::CurveType::Circle)
    {
        Vec3 a = m_points[0], c = m_points[1], b = m_points[2];
        m_centre = c;
        m_a = m_b = (a - c).norm();

        m_e1 = (a - c) / m_a;
        Vec3 w = (b - c) - (b - c).dot(m_e1) * m_e1;
        if (w.norm() < 1e-12 * m_a)
            throw std::invalid_argument("Mesh: circle arc " + std::to_string(curveID) +
                                        " must be shorter than 180 degrees.");
        m_e2 = w.normalized();
        m_t0 = 0.0;
        m_t1 = std::atan2((b - c).dot(m_e2), (b - c).dot(m_e1));
    }
    else if (m_type == Geometry::CurveType::Ellipse)
    {
        Vec3 s = m_points[0], c = m_points[1], major = m_points[2], e = m_points[3];
        m_centre = c;
        m_e1 = (major - c).normalized();

        Vec3 n = (s - c).cross(e - c);
        if (n.norm() < 1e-14)
            n = (s - c).cross(m_e1);
        if (n.norm() < 1e-14)
            n = Vec3(0.0, 0.0, 1.0);
        m_e2 = n.normalized().cross(m_e1).normalized();

        double xs = (s - c).dot(m_e1), ys = (s - c).dot(m_e2);
        double xe = (e - c).dot(m_e1), ye = (e - c).dot(m_e2);

        // x^2 / a^2 + y^2 / b^2 = 1 at the start and end points
        Eigen::Matrix2d A;
        A << xs * xs, ys * ys, xe * xe, ye * ye;
        Eigen::Vector2d inv;

        if (std::abs(A.determinant()) > 1e-14 * A.cwiseAbs().maxCoeff() * A.cwiseAbs().maxCoeff())
            inv = A.partialPivLu().solve(Eigen::Vector2d(1.0, 1.0));
        else
        {
            // Start and end at the same "height": circle through both
            double r = (s - c).norm();
            inv << 1.0 / (r * r), 1.0 / (r * r);
        }

        if ((inv(0) <= 0.0) || (inv(1) <= 0.0))
            throw std::invalid_argument("Mesh: ellipse arc " + std::to_string(curveID) +
                                        " has no ellipse through its start and end points.");

        m_a = 1.0 / std::sqrt(inv(0));
        m_b = 1.0 / std::sqrt(inv(1));
        m_t0 = std::atan2(ys / m_b, xs / m_a);
        m_t1 = std::atan2(ye / m_b, xe / m_a);

        // Shorter arc
        double d = m_t1 - m_t0;
        while (d > pi)
            d -= 2.0 * pi;
        while (d < -pi)
            d += 2.0 * pi;
        m_t1 = m_t0 + d;
    }

    // Arc length table
    m_t.resize(samplesPerCurve + 1);
    m_s.resize(samplesPerCurve + 1);
    Vec3 prev = eval(0.0);
    m_t[0] = 0.0;
    m_s[0] = 0.0;
    for (int i = 1; i <= samplesPerCurve; i++)
    {
        m_t[i] = double(i) / samplesPerCurve;
        Vec3 p = eval(m_t[i]);
        m_s[i] = m_s[i - 1] + (p - prev).norm();
        prev = p;
    }
}

Vec3 CurveShape::evalSpline(double t) const
{
    auto m = static_cast<int>(m_points.size()) - 1;
    if (m == 1)
        return (1.0 - t) * m_points[0] + t * m_points[1];

    int seg = std::min(static_cast<int>(t * m), m - 1);
    double u = t * m - seg;

    auto at = [&](int i) -> Vec3 {
        if (m_closed)
            return m_points[((i % m) + m) % m];
        if (i < 0)
            return 2.0 * m_points[0] - m_points[1];
        if (i > m)
            return 2.0 * m_points[m] - m_points[m - 1];
        return m_points[i];
    };

    Vec3 P0 = at(seg - 1), P1 = at(seg), P2 = at(seg + 1), P3 = at(seg + 2);

    // Catmull-Rom
    return 0.5 * ((2.0 * P1) + (-P0 + P2) * u + (2.0 * P0 - 5.0 * P1 + 4.0 * P2 - P3) * u * u +
                  (-P0 + 3.0 * P1 - 3.0 * P2 + P3) * u * u * u);
}

Vec3 CurveShape::evalBSpline(double t) const
{
    auto m = static_cast<int>(m_points.size()) - 1;
    int p = std::min(3, m);

    // Clamped uniform knot vector
    std::vector<double> knots;
    for (int i = 0; i <= p; i++)
        knots.push_back(0.0);
    for (int i = 1; i <= m - p; i++)
        knots.push_back(double(i) / (m - p + 1));
    for (int i = 0; i <= p; i++)
        knots.push_back(1.0);

    // Knot span
    int k = p;
    while ((k < m) && (t >= knots[k + 1]))
        k++;

    // de Boor
    std::vector<Vec3> d(p + 1);
    for (int j = 0; j <= p; j++)
        d[j] = m_points[j + k - p];

    for (int r = 1; r <= p; r++)
        for (int j = p; j >= r; j--)
        {
            double denom = knots[j + 1 + k - r] - knots[j + k - p];
            double alpha = denom > 0.0 ? (t - knots[j + k - p]) / denom : 0.0;
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }

    return d[p];
}

Vec3 CurveShape::eval(double t) const
{
    switch (m_type)
    {
    case Geometry::CurveType::Spline:
        return evalSpline(t);
    case Geometry::CurveType::BSpline:
        return evalBSpline(t);
    case Geometry::CurveType::Circle:
    case Geometry::CurveType::Ellipse: {
        double angle = m_t0 + t * (m_t1 - m_t0);
        return m_centre + m_a * std::cos(angle) * m_e1 + m_b * std::sin(angle) * m_e2;
    }
    }
    return Vec3::Zero();
}

double CurveShape::length() const
{
    return m_s.back();
}

Vec3 CurveShape::atLength(double s) const
{
    if (s <= 0.0)
        return eval(0.0);
    if (s >= 1.0)
        return eval(1.0);

    double target = s * m_s.back();
    auto it = std::lower_bound(m_s.begin(), m_s.end(), target);
    auto i = static_cast<size_t>(it - m_s.begin());
    double f = (target - m_s[i - 1]) / (m_s[i] - m_s[i - 1]);
    return eval(m_t[i - 1] + f * (m_t[i] - m_t[i - 1]));
}

std::vector<double> curveDivision(const Geometry::Curve &curve, double length, double h0, double h1, bool even)
{
    if (curve.el_on_curve > 0)
    {
        int n = curve.el_on_curve;
        double v = curve.el_distrib_val;

        if ((curve.el_distrib_type == "Progression") && (v > 0.0) && (std::abs(v - 1.0) > 1e-12))
        {
            // Element lengths grow by the factor v
            std::vector<double> xi(n + 1);
            for (int k = 0; k <= n; k++)
                xi[k] = (std::pow(v, k) - 1.0) / (std::pow(v, n) - 1.0);
            return xi;
        }

        if ((curve.el_distrib_type == "Bump") && (v > 0.0))
        {
            // Element size v (relative) at both ends and 1 in the middle (approximates Gmsh's Bump)
            return equalIntegral(n, [v](double x) { return 1.0 / (v + (1.0 - v) * 4.0 * x * (1.0 - x)); });
        }

        std::vector<double> xi(n + 1);
        for (int k = 0; k <= n; k++)
            xi[k] = double(k) / n;
        return xi;
    }

    // From the element sizes at the ends, varying linearly along the curve
    auto h = [h0, h1](double x) { return h0 + (h1 - h0) * x; };

    double integral = 0.0;
    const int m = 200;
    for (int i = 0; i < m; i++)
        integral += 0.5 * (1.0 / h(double(i) / m) + 1.0 / h(double(i + 1) / m)) / m;
    integral *= length;

    int n = std::max(1, static_cast<int>(std::lround(integral)));
    if (even && (n % 2 == 1))
        n++;
    return equalIntegral(n, [&h](double x) { return 1.0 / h(x); });
}

} // namespace calfem::detail
