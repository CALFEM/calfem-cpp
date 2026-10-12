//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// calfem-python drawing functions for 2D elements (eldraw2, eldisp2, dispbeam2, secforce2,
// scalgraph2), ported from calfem/vis_mpl.py.

#include <calfem/vis.h>

#include <calfem/elements/beam.h>

#include "detail/vis_codes.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <tuple>

namespace calfem::vis {

namespace {

// ---- plotpar codes ----

using codes::colorCode;
using codes::lineCode;

// Node marks as in calfem-python's pltstyle2 (eldisp2, dispbeam2)
Marker nodeCode(int code, const std::string &what)
{
    switch (code)
    {
    case 0:
        return Marker::None;
    case 1:
        return Marker::Circle;
    case 2:
        return Marker::Star;
    case 3:
        return Marker::Dot;
    default:
        throw std::invalid_argument("Invalid value for " + what + " (node mark 0-3), got " + std::to_string(code) + ".");
    }
}

// ---- Helpers ----

void checkElementCoords(const Matrix &ex, const Matrix &ey, const char *caller)
{
    if (ex.rows() != ey.rows() || ex.cols() != ey.cols())
        throw std::invalid_argument(std::string(caller) + ": check size of ex, ey dimensions.");
}

// Points of row i of x, y in the given column order
std::vector<Point> rowPoints(const Matrix &x, const Matrix &y, Eigen::Index i, const std::vector<int> &order)
{
    std::vector<Point> points;
    points.reserve(order.size());
    for (int j : order)
        points.push_back({x(i, j), y(i, j)});
    return points;
}

std::vector<int> range(Eigen::Index n)
{
    std::vector<int> r(static_cast<std::size_t>(n));
    for (Eigen::Index i = 0; i < n; i++)
        r[static_cast<std::size_t>(i)] = int(i);
    return r;
}

// One line per element (calfem-python's draw_elements)
void addElementLines(Figure &fig, const Matrix &x, const Matrix &y, const std::vector<int> &order, const Color &color,
                     LineStyle style, bool closed)
{
    for (Eigen::Index i = 0; i < x.rows(); i++)
    {
        Line line;
        line.points = rowPoints(x, y, i, order);
        line.color = color;
        line.style = style;
        line.width = 1.0;
        line.closed = closed;
        fig.add(std::move(line));
    }
}

// Markers at all element nodes (calfem-python's draw_node_circles)
void addNodeMarkers(Figure &fig, const Matrix &x, const Matrix &y, Marker marker, const Color &color)
{
    if (marker == Marker::None)
        return;
    Markers markers;
    markers.marker = marker;
    markers.color = color;
    for (Eigen::Index i = 0; i < x.rows(); i++)
        for (Eigen::Index j = 0; j < x.cols(); j++)
            markers.points.push_back({x(i, j), y(i, j)});
    fig.add(std::move(markers));
}

void addPolyline(Figure &fig, std::vector<Point> points, const Color &color, double width,
                 LineStyle style = LineStyle::Solid)
{
    Line line;
    line.points = std::move(points);
    line.color = color;
    line.width = width;
    line.style = style;
    fig.add(std::move(line));
}

// numpy.linspace(0, L, n)
std::vector<double> linspace(double L, Eigen::Index n)
{
    std::vector<double> v(static_cast<std::size_t>(n));
    if (n == 1)
        v[0] = 0.0;
    else
    {
        double step = L / double(n - 1);
        for (Eigen::Index i = 0; i < n; i++)
            v[static_cast<std::size_t>(i)] = double(i) * step;
        v.back() = L;
    }
    return v;
}

// Python's str() of a float (shortest representation that reads back the same value)
std::string pythonFloat(double value)
{
    if (std::isnan(value))
        return "nan";
    if (std::isinf(value))
        return value > 0 ? "inf" : "-inf";
    if (value == 0.0)
        return std::signbit(value) ? "-0.0" : "0.0";

    char buffer[64];
    int precision = 1;
    for (; precision <= 17; precision++)
    {
        std::snprintf(buffer, sizeof(buffer), "%.*e", precision - 1, value);
        if (std::strtod(buffer, nullptr) == value)
            break;
    }

    // buffer: [-]d.ddde[+-]xx
    std::string s(buffer);
    bool negative = s[0] == '-';
    if (negative)
        s.erase(0, 1);
    auto ePos = s.find('e');
    int exponent = std::atoi(s.c_str() + ePos + 1);
    std::string digits = s.substr(0, ePos);
    digits.erase(std::remove(digits.begin(), digits.end(), '.'), digits.end());
    while (digits.size() > 1 && digits.back() == '0')
        digits.pop_back();

    std::string result;
    if (exponent >= -4 && exponent < 16)
    {
        if (exponent < 0)
            result = "0." + std::string(std::size_t(-exponent - 1), '0') + digits;
        else if (std::size_t(exponent) + 1 >= digits.size())
            result = digits + std::string(std::size_t(exponent) + 1 - digits.size(), '0') + ".0";
        else
            result = digits.substr(0, std::size_t(exponent) + 1) + "." + digits.substr(std::size_t(exponent) + 1);
    }
    else
    {
        result = digits.substr(0, 1);
        if (digits.size() > 1)
            result += "." + digits.substr(1);
        std::snprintf(buffer, sizeof(buffer), "e%c%02d", exponent < 0 ? '-' : '+', std::abs(exponent));
        result += buffer;
    }
    return negative ? "-" + result : result;
}

// Element direction (unit vector) and length of a 2-node element
void elementAxis(const RowVec &ex, const RowVec &ey, double &L, double &nx, double &ny, const char *caller)
{
    if (ex.size() != 2 || ey.size() != 2)
        throw std::invalid_argument(std::string(caller) + ": ex and ey must hold the 2 node coordinates of one element.");
    double dx = ex(1) - ex(0), dy = ey(1) - ey(0);
    L = std::sqrt(dx * dx + dy * dy);
    nx = dx / L;
    ny = dy / L;
}

// 8-node element boundary through n points per edge (calfem-python's _quad8_edges)
std::tuple<Matrix, Matrix> quad8Edges(const Matrix &x, const Matrix &y, int n = 7)
{
    auto shape = [](double t, double s) {
        RowVec N(8);
        N << -0.25 * (1 - t) * (1 - s) * (1 + t + s), -0.25 * (1 + t) * (1 - s) * (1 - t + s),
            -0.25 * (1 + t) * (1 + s) * (1 - t - s), -0.25 * (1 - t) * (1 + s) * (1 + t - s),
            0.5 * (1 - t * t) * (1 - s), 0.5 * (1 + t) * (1 - s * s), 0.5 * (1 - t * t) * (1 + s),
            0.5 * (1 - t) * (1 - s * s);
        return N;
    };

    std::vector<double> g(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++)
        g[std::size_t(i)] = -1.0 + 2.0 * i / double(n - 1);

    std::vector<std::pair<double, double>> path;
    for (int i = 0; i < n - 1; i++)
        path.push_back({g[std::size_t(i)], -1.0});
    for (int i = 0; i < n - 1; i++)
        path.push_back({1.0, g[std::size_t(i)]});
    for (int i = n - 1; i > 0; i--)
        path.push_back({g[std::size_t(i)], 1.0});
    for (int i = n - 1; i >= 0; i--)
        path.push_back({-1.0, g[std::size_t(i)]});

    Matrix N(Eigen::Index(path.size()), 8);
    for (std::size_t k = 0; k < path.size(); k++)
        N.row(Eigen::Index(k)) = shape(path[k].first, path[k].second);

    return {x * N.transpose(), y * N.transpose()};
}

} // namespace

// ---- eldraw2 ----

void eldraw2(const Matrix &ex, const Matrix &ey, const std::array<int, 3> &plotpar, const IntRowVec &elnum)
{
    eldraw2(gcf(), ex, ey, plotpar, elnum);
}

void eldraw2(Figure &fig, const Matrix &ex, const Matrix &ey, const std::array<int, 3> &plotpar,
             const IntRowVec &elnum)
{
    checkElementCoords(ex, ey, "eldraw2");

    LineStyle style = lineCode(plotpar[0], "plotpar[0]");
    Color color = colorCode(plotpar[1], "plotpar[1]");

    Marker marker;
    switch (plotpar[2])
    {
    case 0:
        marker = Marker::None;
        break;
    case 1:
        marker = Marker::Circle;
        break;
    case 2:
        marker = Marker::Cross;
        break;
    default:
        throw std::invalid_argument("eldraw2: invalid value for plotpar[2] (node mark 0-2), got " +
                                    std::to_string(plotpar[2]) + ".");
    }

    if (elnum.size() > 0 && elnum.size() != ex.rows())
        throw std::invalid_argument("eldraw2: elnum must contain one number per element.");

    auto nen = ex.cols();

    // 8-node elements: corner and midside nodes in order along the boundary
    std::vector<int> order = nen == 8 ? std::vector<int>{0, 4, 1, 5, 2, 6, 3, 7} : range(nen);

    addElementLines(fig, ex, ey, order, color, style, nen > 2);
    addNodeMarkers(fig, ex, ey, marker, color);

    for (Eigen::Index i = 0; i < elnum.size(); i++)
    {
        Text t;
        t.position = {ex.row(i).mean(), ey.row(i).mean()};
        t.text = std::to_string(elnum(i));
        t.color = color;
        t.halign = HAlign::Center;
        t.valign = VAlign::Center;
        fig.add(std::move(t));
    }
}

// ---- eldisp2 ----

double eldisp2(const Matrix &ex, const Matrix &ey, const Matrix &ed, const std::array<int, 3> &plotpar,
               std::optional<double> sfac)
{
    return eldisp2(gcf(), ex, ey, ed, plotpar, sfac);
}

double eldisp2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed, const std::array<int, 3> &plotpar,
               std::optional<double> sfac)
{
    checkElementCoords(ex, ey, "eldisp2");
    if (ed.rows() != ex.rows())
        throw std::invalid_argument("eldisp2: check size of ed/ex dimensions.");

    LineStyle style = lineCode(plotpar[0], "plotpar[0]");
    Color color = colorCode(plotpar[1], "plotpar[1]");
    Marker marker = nodeCode(plotpar[2], "plotpar[2]");

    double dlMax = std::max(ex.maxCoeff() - ex.minCoeff(), ey.maxCoeff() - ey.minCoeff());
    double edMax = ed.cwiseAbs().maxCoeff();
    double k = sfac ? *sfac : (edMax > 0.0 ? 0.1 * dlMax / edMax : 1.0);

    auto nen = ex.cols(), ned = ed.cols();
    Matrix x, y, xc, yc;
    bool closed = false;

    // Displacement components at the nodes (every other column of ed)
    auto nodeDisp = [&ed](Eigen::Index offset, Eigen::Index stride, Eigen::Index n) {
        Matrix d(ed.rows(), n);
        for (Eigen::Index j = 0; j < n; j++)
            d.col(j) = ed.col(offset + j * stride);
        return d;
    };

    if (nen == 2 && ned == 4) // bars
    {
        x = ex + k * nodeDisp(0, 2, 2);
        y = ey + k * nodeDisp(1, 2, 2);
        xc = x;
        yc = y;
    }
    else if (nen == 2 && ned == 6) // beams
    {
        x = ex + k * nodeDisp(0, 3, 2);
        y = ey + k * nodeDisp(1, 3, 2);
        std::tie(xc, yc) = beam2crd(ex, ey, ed, k);
    }
    else if ((nen == 3 || nen == 4) && ned == 2 * nen) // triangles and quadrilaterals
    {
        x = ex + k * nodeDisp(0, 2, nen);
        y = ey + k * nodeDisp(1, 2, nen);
        xc = x;
        yc = y;
        closed = true;
    }
    else if (nen == 8 && ned == 16) // 8-node isoparametric elements
    {
        x = ex + k * nodeDisp(0, 2, 8);
        y = ey + k * nodeDisp(1, 2, 8);
        std::tie(xc, yc) = quad8Edges(x, y);
    }
    else
        throw std::invalid_argument("eldisp2: element type is not supported (ex: " + std::to_string(nen) +
                                    " columns, ed: " + std::to_string(ned) + " columns).");

    addElementLines(fig, xc, yc, range(xc.cols()), color, style, closed);
    addNodeMarkers(fig, x, y, marker, Color{});

    return k;
}

// ---- dispbeam2 ----

double dispbeam2(const RowVec &ex, const RowVec &ey, const Matrix &edi, const std::array<int, 3> &plotpar,
                 std::optional<double> sfac)
{
    return dispbeam2(gcf(), ex, ey, edi, plotpar, sfac);
}

double dispbeam2(Figure &fig, const RowVec &ex, const RowVec &ey, const Matrix &edi,
                 const std::array<int, 3> &plotpar, std::optional<double> sfac)
{
    double L, nx, ny;
    elementAxis(ex, ey, L, nx, ny, "dispbeam2");
    if (edi.cols() != 2 || edi.rows() < 1)
        throw std::invalid_argument("dispbeam2: check size of edi dimension (one row [u, v] per point).");

    LineStyle style = lineCode(plotpar[0], "plotpar[0]");
    Color color = colorCode(plotpar[1], "plotpar[1]");
    Marker marker = nodeCode(plotpar[2], "plotpar[2]");

    double ediMax = edi.cwiseAbs().maxCoeff();
    double k = sfac ? *sfac : (ediMax > 0.0 ? 0.1 * L / ediMax : 1.0);

    auto n = edi.rows();
    auto eci = linspace(L, n);

    std::vector<Point> points;
    for (Eigen::Index i = 0; i < n; i++)
    {
        double x = ex(0) + eci[std::size_t(i)] * nx, y = ey(0) + eci[std::size_t(i)] * ny;
        double u = k * edi(i, 0), v = k * edi(i, 1);
        points.push_back({x + u * nx - v * ny, y + u * ny + v * nx});
    }

    Markers markers;
    markers.marker = marker;
    markers.points = {points.front(), points.back()};

    addPolyline(fig, std::move(points), color, 1.0, style);
    if (marker != Marker::None)
        fig.add(std::move(markers));

    return k;
}

// ---- secforce2 ----

double secforce2(const RowVec &ex, const RowVec &ey, const Matrix &es, const std::array<int, 2> &plotpar,
                 std::optional<double> sfac, const Matrix &eci)
{
    return secforce2(gcf(), ex, ey, es, plotpar, sfac, eci);
}

double secforce2(Figure &fig, const RowVec &ex, const RowVec &ey, const Matrix &es, const std::array<int, 2> &plotpar,
                 std::optional<double> sfac, const Matrix &eci)
{
    double L, nx, ny;
    elementAxis(ex, ey, L, nx, ny, "secforce2");
    if ((es.rows() != 1 && es.cols() != 1) || es.size() == 0)
        throw std::invalid_argument("secforce2: es must be a vector with one section force per evaluation point.");

    Color color = colorCode(plotpar[0], "plotpar[0]");
    Color elementColor = colorCode(plotpar[1], "plotpar[1]");

    auto n = es.size();
    std::vector<double> xi;
    if (eci.size() == 0)
        xi = linspace(L, n);
    else if (eci.size() != n || (eci.rows() != 1 && eci.cols() != 1))
        throw std::invalid_argument("secforce2: check size of eci dimension (one coordinate per value of es).");
    else
        xi.assign(eci.data(), eci.data() + n);

    double esMax = es.cwiseAbs().maxCoeff();
    double k = sfac ? *sfac : (esMax > 0.0 ? 0.2 * L / esMax : 1.0);

    std::vector<Point> base, diagram;
    for (Eigen::Index i = 0; i < n; i++)
    {
        Point b{ex(0) + xi[std::size_t(i)] * nx, ey(0) + xi[std::size_t(i)] * ny};
        double s = k * es(i);
        base.push_back(b);
        diagram.push_back({b.x + s * ny, b.y - s * nx});
    }

    addPolyline(fig, diagram, color, 1.0);
    for (Eigen::Index i = 0; i < n; i++)
        addPolyline(fig, {base[std::size_t(i)], diagram[std::size_t(i)]}, color, 1.0);
    addPolyline(fig, {{ex(0), ey(0)}, {ex(1), ey(1)}}, elementColor, 2.0);

    return k;
}

// ---- scalgraph2 ----

void scalgraph2(double sfac, const std::vector<double> &magnitude, int plotpar)
{
    scalgraph2(gcf(), sfac, magnitude, plotpar);
}

void scalgraph2(Figure &fig, double sfac, const std::vector<double> &magnitude, int plotpar)
{
    if (magnitude.size() != 1 && magnitude.size() != 3)
        throw std::invalid_argument("scalgraph2: magnitude must be {N} or {N, x, y}.");

    double N = magnitude[0];
    double x = magnitude.size() == 3 ? magnitude[1] : 0.0;
    double y = magnitude.size() == 3 ? magnitude[2] : -0.5;
    double L = N * sfac;

    Color color = colorCode(plotpar, "plotpar");

    addPolyline(fig, {{x, y}, {x + L, y}}, color, 1.0);
    addPolyline(fig, {{x, y - L / 20.0}, {x, y + L / 20.0}}, color, 1.0);
    addPolyline(fig, {{x + L, y - L / 20.0}, {x + L, y + L / 20.0}}, color, 1.0);

    Text t;
    t.position = {x + L * 1.1, y - L / 20.0};
    t.text = pythonFloat(N);
    fig.add(std::move(t));
}

} // namespace calfem::vis
