//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Geometry description, ported from calfem-python (geometry.py).

#include <calfem/geometry.h>

#include <algorithm>
#include <cctype>
#include <limits>
#include <stdexcept>

namespace calfem {

namespace {

std::string titleCase(const std::string &s)
{
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (!t.empty())
        t[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(t[0])));
    return t;
}

template <typename T> T &lookup(std::map<int, T> &m, int ID, const char *what)
{
    auto it = m.find(ID);
    if (it == m.end())
        throw std::invalid_argument(std::string("Geometry: ") + what + " " + std::to_string(ID) + " does not exist.");
    return it->second;
}

} // namespace

template <typename T> int Geometry::smallestFreeKey(const std::map<int, T> &m)
{
    int i = 0;
    for (const auto &[key, value] : m)
    {
        if (key != i)
            return i;
        i++;
    }
    return i;
}

// ---- Points ----

void Geometry::point(const std::vector<double> &coord, int ID, int marker, double el_size)
{
    if ((coord.size() != 2) && (coord.size() != 3))
        throw std::invalid_argument("Geometry::point: coord must have 2 or 3 components.");

    if (coord.size() == 3)
        m_is3D = true;

    if (ID < 0)
        ID = m_pointIDspecified ? smallestFreeKey(m_points) : m_nextPointID++;
    else
        m_pointIDspecified = true;

    Point p;
    p.coord << coord[0], coord[1], coord.size() == 3 ? coord[2] : 0.0;
    p.el_size = el_size;
    p.marker = marker;
    m_points[ID] = p;
}

// ---- Curves ----

void Geometry::addCurve(CurveType type, const std::vector<int> &points, int ID, int marker, int el_on_curve,
                        const std::string &el_distrib_type, double el_distrib_val)
{
    if (ID < 0)
        ID = m_curveIDspecified ? smallestFreeKey(m_curves) : m_nextCurveID++;
    else
        m_curveIDspecified = true;

    Curve c;
    c.type = type;
    c.points = points;
    c.marker = marker;
    c.el_on_curve = el_on_curve;
    c.el_distrib_val = el_distrib_val;

    if (!el_distrib_type.empty())
    {
        c.el_distrib_type = titleCase(el_distrib_type);
        if ((c.el_distrib_type != "Bump") && (c.el_distrib_type != "Progression"))
            throw std::invalid_argument("Geometry: el_distrib_type must be \"bump\" or \"progression\" (curve " +
                                        std::to_string(ID) + ").");
    }

    m_curves[ID] = c;
}

void Geometry::spline(const std::vector<int> &points, int ID, int marker, int el_on_curve,
                      const std::string &el_distrib_type, double el_distrib_val)
{
    if (points.size() < 2)
        throw std::invalid_argument("Geometry::spline: at least two points are needed.");
    addCurve(CurveType::Spline, points, ID, marker, el_on_curve, el_distrib_type, el_distrib_val);
}

void Geometry::bspline(const std::vector<int> &points, int ID, int marker, int el_on_curve,
                       const std::string &el_distrib_type, double el_distrib_val)
{
    if (points.size() < 2)
        throw std::invalid_argument("Geometry::bspline: at least two points are needed.");
    addCurve(CurveType::BSpline, points, ID, marker, el_on_curve, el_distrib_type, el_distrib_val);
}

void Geometry::circle(const std::vector<int> &points, int ID, int marker, int el_on_curve,
                      const std::string &el_distrib_type, double el_distrib_val)
{
    if (points.size() != 3)
        throw std::invalid_argument("Geometry::circle: points must be [start, centre, end].");
    addCurve(CurveType::Circle, points, ID, marker, el_on_curve, el_distrib_type, el_distrib_val);
}

void Geometry::ellipse(const std::vector<int> &points, int ID, int marker, int el_on_curve,
                       const std::string &el_distrib_type, double el_distrib_val)
{
    if (points.size() != 4)
        throw std::invalid_argument("Geometry::ellipse: points must be [start, centre, major axis point, end].");
    addCurve(CurveType::Ellipse, points, ID, marker, el_on_curve, el_distrib_type, el_distrib_val);
}

// ---- Surfaces ----

void Geometry::addSurface(const std::string &name, const std::vector<int> &outer_loop,
                          const std::vector<std::vector<int>> &holes, int ID, int marker, bool is_structured)
{
    if (ID < 0)
        ID = m_surfaceIDspecified ? smallestFreeKey(m_surfaces) : m_nextSurfaceID++;
    else
        m_surfaceIDspecified = true;

    Surface s;
    s.name = name;
    s.outer_loop = outer_loop;
    s.holes = holes;
    s.marker = marker;
    s.is_structured = is_structured;
    m_surfaces[ID] = s;
}

void Geometry::surface(const std::vector<int> &outer_loop, const std::vector<std::vector<int>> &holes, int ID,
                       int marker)
{
    addSurface("Plane Surface", outer_loop, holes, ID, marker, false);
}

void Geometry::ruled_surface(const std::vector<int> &outer_loop, int ID, int marker)
{
    if ((outer_loop.size() != 3) && (outer_loop.size() != 4))
        throw std::invalid_argument("Geometry::ruled_surface: outer_loop must have 3 or 4 curves.");
    addSurface("Surface", outer_loop, {}, ID, marker, false);
}

void Geometry::struct_surface(const std::vector<int> &outer_loop, int ID, int marker)
{
    if (outer_loop.size() != 4)
        throw std::invalid_argument("Geometry::struct_surface: outer_loop must have 4 curves.");

    int n[4];
    for (int i = 0; i < 4; i++)
    {
        auto it = m_curves.find(outer_loop[i]);
        if (it == m_curves.end())
            throw std::invalid_argument("Geometry::struct_surface: curve " + std::to_string(outer_loop[i]) +
                                        " does not exist.");
        n[i] = it->second.el_on_curve;
        if (n[i] < 1)
            throw std::invalid_argument("Geometry::struct_surface: curve " + std::to_string(outer_loop[i]) +
                                        " has no el_on_curve (needed for structured surfaces).");
    }

    if ((n[0] != n[2]) || (n[1] != n[3]))
        throw std::invalid_argument(
            "Geometry::struct_surface: opposite curves must have the same number of elements (el_on_curve).");

    addSurface("Surface", outer_loop, {}, ID, marker, true);
}

// ---- Volumes ----

void Geometry::addVolume(const std::vector<int> &outer_surfaces, const std::vector<std::vector<int>> &holes, int ID,
                         int marker, bool is_structured)
{
    if (ID < 0)
        ID = m_volumeIDspecified ? smallestFreeKey(m_volumes) : m_nextVolumeID++;
    else
        m_volumeIDspecified = true;

    Volume v;
    v.outer_surfaces = outer_surfaces;
    v.holes = holes;
    v.marker = marker;
    v.is_structured = is_structured;
    m_volumes[ID] = v;
}

void Geometry::volume(const std::vector<int> &outer_surfaces, const std::vector<std::vector<int>> &holes, int ID,
                      int marker)
{
    addVolume(outer_surfaces, holes, ID, marker, false);
}

void Geometry::struct_volume(const std::vector<int> &outer_surfaces, int ID, int marker)
{
    addVolume(outer_surfaces, {}, ID, marker, true);
}

// ---- Markers ----

void Geometry::setPointMarker(int ID, int marker)
{
    lookup(m_points, ID, "point").marker = marker;
}

void Geometry::setCurveMarker(int ID, int marker)
{
    lookup(m_curves, ID, "curve").marker = marker;
}

void Geometry::setSurfaceMarker(int ID, int marker)
{
    lookup(m_surfaces, ID, "surface").marker = marker;
}

void Geometry::setVolumeMarker(int ID, int marker)
{
    lookup(m_volumes, ID, "volume").marker = marker;
}

void Geometry::pointMarker(int ID, int marker)
{
    setPointMarker(ID, marker);
}

void Geometry::curveMarker(int ID, int marker)
{
    setCurveMarker(ID, marker);
}

void Geometry::surfaceMarker(int ID, int marker)
{
    setSurfaceMarker(ID, marker);
}

// ---- Removal ----

void Geometry::removePoint(int ID)
{
    m_points.erase(ID);
}

void Geometry::removeCurve(int ID)
{
    m_curves.erase(ID);
}

void Geometry::removeSurface(int ID)
{
    m_surfaces.erase(ID);
}

void Geometry::removeVolume(int ID)
{
    m_volumes.erase(ID);
}

// ---- Queries ----

const std::map<int, Geometry::Point> &Geometry::points() const
{
    return m_points;
}

const std::map<int, Geometry::Curve> &Geometry::curves() const
{
    return m_curves;
}

const std::map<int, Geometry::Surface> &Geometry::surfaces() const
{
    return m_surfaces;
}

const std::map<int, Geometry::Volume> &Geometry::volumes() const
{
    return m_volumes;
}

bool Geometry::is3D() const
{
    return m_is3D;
}

Matrix Geometry::getPointCoords() const
{
    Matrix c(static_cast<Eigen::Index>(m_points.size()), 3);
    Eigen::Index i = 0;
    for (const auto &[ID, p] : m_points)
        c.row(i++) = p.coord;
    return c;
}

Matrix Geometry::getPointCoords(const std::vector<int> &IDs) const
{
    Matrix c(static_cast<Eigen::Index>(IDs.size()), 3);
    for (size_t i = 0; i < IDs.size(); i++)
    {
        auto it = m_points.find(IDs[i]);
        if (it == m_points.end())
            throw std::invalid_argument("Geometry: point " + std::to_string(IDs[i]) + " does not exist.");
        c.row(static_cast<Eigen::Index>(i)) = it->second.coord;
    }
    return c;
}

std::tuple<double, double, double, double> Geometry::bounding_box_2d() const
{
    double min_x = 1e300, max_x = -1e300, min_y = 1e300, max_y = -1e300;
    for (const auto &[ID, p] : m_points)
    {
        min_x = std::min(min_x, p.coord(0));
        max_x = std::max(max_x, p.coord(0));
        min_y = std::min(min_y, p.coord(1));
        max_y = std::max(max_y, p.coord(1));
    }
    return {min_x, max_x, min_y, max_y};
}

std::vector<int> Geometry::pointsOnCurves(const std::vector<int> &IDs) const
{
    std::set<int> points;
    for (int ID : IDs)
    {
        auto it = m_curves.find(ID);
        if (it != m_curves.end())
            points.insert(it->second.points.begin(), it->second.points.end());
    }
    return {points.begin(), points.end()};
}

std::tuple<std::vector<int>, std::vector<int>> Geometry::stuffOnSurfaces(const std::vector<int> &IDs) const
{
    std::set<int> curves;
    for (int ID : IDs)
    {
        auto it = m_surfaces.find(ID);
        if (it == m_surfaces.end())
            continue;
        curves.insert(it->second.outer_loop.begin(), it->second.outer_loop.end());
        for (const auto &hole : it->second.holes)
            curves.insert(hole.begin(), hole.end());
    }
    std::vector<int> curveList(curves.begin(), curves.end());
    return {pointsOnCurves(curveList), curveList};
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> Geometry::stuffOnVolumes(
    const std::vector<int> &IDs) const
{
    std::set<int> surfaces;
    for (int ID : IDs)
    {
        auto it = m_volumes.find(ID);
        if (it == m_volumes.end())
            continue;
        surfaces.insert(it->second.outer_surfaces.begin(), it->second.outer_surfaces.end());
        for (const auto &hole : it->second.holes)
            surfaces.insert(hole.begin(), hole.end());
    }
    std::vector<int> surfaceList(surfaces.begin(), surfaces.end());
    auto [pointList, curveList] = stuffOnSurfaces(surfaceList);
    return {pointList, curveList, surfaceList};
}

} // namespace calfem
