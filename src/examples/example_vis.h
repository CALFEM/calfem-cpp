//
// CALFEM C++ examples: helpers shared by the *_vis examples
//
// Graphs (x-y curves), arrows and projection of 3D models built from calfem::vis items, and the
// common ending of the examples: show the figures in Qt windows when the Qt backend is built
// (CALFEM_EXAMPLE_QT), otherwise or with the argument --svg write them to <prefix>_<name>.svg.
//
// SPDX-License-Identifier: MIT
//

#pragma once

#include <calfem/calfem.h>

#ifdef CALFEM_EXAMPLE_QT
#include <calfem/vis_qt.h>
#endif

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace example {

namespace vis = calfem::vis;

const vis::Color colorBlack{0.0, 0.0, 0.0};
const vis::Color colorBlue{0.0, 0.0, 1.0};
const vis::Color colorRed{1.0, 0.0, 0.0};
const vis::Color colorGreen{0.0, 0.6, 0.0};
const vis::Color colorGray{0.6, 0.6, 0.6};

// Makes fig a graph: independent x and y scaling and a grid.
inline void graph(vis::Figure &fig)
{
    fig.set_equal_aspect(false);
    fig.set_grid(true);
}

// Curve through (x[i], y[i]).
inline void plot(vis::Figure &fig, const std::vector<double> &x, const std::vector<double> &y,
                 const vis::Color &color = colorBlue, vis::LineStyle style = vis::LineStyle::Solid, double width = 1.5)
{
    vis::Line line;
    for (std::size_t i = 0; i < x.size() && i < y.size(); i++)
        line.points.push_back({x[i], y[i]});
    line.color = color;
    line.style = style;
    line.width = width;
    fig.add(line);
}

// Markers at (x[i], y[i]).
inline void plotMarkers(vis::Figure &fig, const std::vector<double> &x, const std::vector<double> &y,
                        const vis::Color &color = colorRed, vis::Marker marker = vis::Marker::Circle, double size = 6.0)
{
    vis::Markers markers;
    for (std::size_t i = 0; i < x.size() && i < y.size(); i++)
        markers.points.push_back({x[i], y[i]});
    markers.marker = marker;
    markers.color = color;
    markers.size = size;
    fig.add(markers);
}

// Arrow from (x, y) to (x + dx, y + dy); the head is 25% of the arrow length.
inline void arrow(vis::Figure &fig, double x, double y, double dx, double dy, const vis::Color &color = colorBlue)
{
    double len = std::hypot(dx, dy);
    if (len == 0.0)
        return;
    double ux = dx / len, uy = dy / len, h = 0.25 * len;
    vis::Line shaft{{{x, y}, {x + dx, y + dy}}, color};
    vis::Line head{{{x + dx - h * (ux - 0.5 * uy), y + dy - h * (uy + 0.5 * ux)},
                    {x + dx, y + dy},
                    {x + dx - h * (ux + 0.5 * uy), y + dy - h * (uy - 0.5 * ux)}},
                   color};
    fig.add(shaft);
    fig.add(head);
}

// Orthographic projection of (x, y, z), z pointing up, seen from azimuth 25 degrees (from the
// x-axis towards the y-axis) and elevation 25 degrees. (The isometric view, 45 degrees azimuth,
// would hide bars along the line x = y behind each other.)
inline vis::Point project(double x, double y, double z)
{
    const double deg = 3.14159265358979323846 / 180.0;
    const double ca = std::cos(25.0 * deg), sa = std::sin(25.0 * deg);
    const double ce = std::cos(25.0 * deg), se = std::sin(25.0 * deg);
    return {-x * sa + y * ca, z * ce - (x * ca + y * sa) * se};
}

// Coordinate axes x, y, z of length len at (x0, y0, z0), projected with project().
inline void axes3d(vis::Figure &fig, double x0, double y0, double z0, double len)
{
    const char *names[] = {"x", "y", "z"};
    for (int d = 0; d < 3; d++)
    {
        double e[3] = {0.0, 0.0, 0.0};
        e[d] = len;
        vis::Point p0 = project(x0, y0, z0), p1 = project(x0 + e[0], y0 + e[1], z0 + e[2]);
        arrow(fig, p0.x, p0.y, p1.x - p0.x, p1.y - p0.y, colorGray);
        vis::text(fig, names[d], {p1.x + 0.05 * len, p1.y + 0.05 * len}, colorGray, 9.0);
    }
}

// Shows all figures (Qt) or writes figure n to <prefix>_<names[n - 1]>.svg.
inline int finish(int argc, char *argv[], const std::string &prefix, const std::vector<std::string> &names)
{
    bool svg = false;
    for (int i = 1; i < argc; i++)
        svg = svg || std::strcmp(argv[i], "--svg") == 0;

#ifdef CALFEM_EXAMPLE_QT
    if (!svg)
    {
        vis::show();
        return 0;
    }
#endif
    (void)svg;

    for (std::size_t n = 1; n <= names.size(); n++)
    {
        std::string filename = prefix + "_" + names[n - 1] + ".svg";
        vis::save_figure(filename, vis::figure(int(n)));
        std::cout << "Figure " << n << " written to " << filename << "\n";
    }
    return 0;
}

} // namespace example
