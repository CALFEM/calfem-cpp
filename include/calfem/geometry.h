//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Geometry description for mesh generation (port of calfem-python's calfem.geometry).

#pragma once

#include <calfem/types.h>

#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace calfem {

// Geometry built from points, curves, surfaces and volumes, as in calfem-python:
//
//   Geometry g;
//   g.point({0.0, 0.0});            // point 0
//   g.point({1.0, 0.0});            // point 1
//   ...
//   g.spline({0, 1}, -1, 10);       // curve 0 from point 0 to 1, marker 10
//   ...
//   g.surface({0, 1, 2, 3});        // surface bounded by curves 0..3
//
// IDs are zero-based. When no ID is given, the next free ID is used (as calfem-python: the
// next number in sequence until an ID is given explicitly, then the smallest free ID).
// Markers identify boundaries and regions in the mesh (bdofs, element_markers).
//
// Curves are defined by point IDs:
//   spline   two or more points; a straight line for two points
//   bspline  control points of a B-spline (passes through the first and last point)
//   circle   [start, centre, end], an arc of less than 180 degrees
//   ellipse  [start, centre, point on the major axis, end]
//
// Element distribution on a curve: el_on_curve elements (-1: from the element sizes), optionally
// graded with el_distrib_type "progression" (element sizes grow by el_distrib_val) or "bump"
// (smaller elements at both ends for el_distrib_val < 1).

class Geometry {
public:
    struct Point {
        RowVec3 coord;
        double el_size = 1.0;
        int marker = 0;
    };

    enum class CurveType { Spline, BSpline, Circle, Ellipse };

    struct Curve {
        CurveType type = CurveType::Spline;
        std::vector<int> points;
        int marker = 0;
        int el_on_curve = -1;        // -1: not specified
        std::string el_distrib_type; // "", "Progression" or "Bump"
        double el_distrib_val = 0.0;
    };

    struct Surface {
        std::string name; // "Plane Surface" or "Surface" (ruled/structured), as in the .geo format
        std::vector<int> outer_loop;
        std::vector<std::vector<int>> holes;
        int marker = 0;
        bool is_structured = false;
    };

    struct Volume {
        std::vector<int> outer_surfaces;
        std::vector<std::vector<int>> holes;
        int marker = 0;
        bool is_structured = false;
    };

    Geometry() = default;

    // ---- Points ----

    // coord with 2 or 3 components; a 3D point makes the geometry three-dimensional
    void point(const std::vector<double> &coord, int ID = -1, int marker = 0, double el_size = 1.0);

    // ---- Curves ----

    void spline(const std::vector<int> &points, int ID = -1, int marker = 0, int el_on_curve = -1,
                const std::string &el_distrib_type = "", double el_distrib_val = 0.0);
    void bspline(const std::vector<int> &points, int ID = -1, int marker = 0, int el_on_curve = -1,
                 const std::string &el_distrib_type = "", double el_distrib_val = 0.0);
    void circle(const std::vector<int> &points, int ID = -1, int marker = 0, int el_on_curve = -1,
                const std::string &el_distrib_type = "", double el_distrib_val = 0.0);
    void ellipse(const std::vector<int> &points, int ID = -1, int marker = 0, int el_on_curve = -1,
                 const std::string &el_distrib_type = "", double el_distrib_val = 0.0);

    // ---- Surfaces ----

    void surface(const std::vector<int> &outer_loop, const std::vector<std::vector<int>> &holes = {}, int ID = -1,
                 int marker = 0);
    void ruled_surface(const std::vector<int> &outer_loop, int ID = -1, int marker = 0);

    // Structured (transfinite) surface bounded by four curves with el_on_curve set, opposite
    // curves with the same number of elements
    void struct_surface(const std::vector<int> &outer_loop, int ID = -1, int marker = 0);

    // ---- Volumes ----

    void volume(const std::vector<int> &outer_surfaces, const std::vector<std::vector<int>> &holes = {}, int ID = -1,
                int marker = 0);
    void struct_volume(const std::vector<int> &outer_surfaces, int ID = -1, int marker = 0);

    // ---- Markers ----

    void setPointMarker(int ID, int marker);
    void setCurveMarker(int ID, int marker);
    void setSurfaceMarker(int ID, int marker);
    void setVolumeMarker(int ID, int marker);

    void pointMarker(int ID, int marker);
    void curveMarker(int ID, int marker);
    void surfaceMarker(int ID, int marker);

    // ---- Removal ----

    void removePoint(int ID);
    void removeCurve(int ID);
    void removeSurface(int ID);
    void removeVolume(int ID);

    // ---- Queries ----

    const std::map<int, Point> &points() const;
    const std::map<int, Curve> &curves() const;
    const std::map<int, Surface> &surfaces() const;
    const std::map<int, Volume> &volumes() const;

    bool is3D() const;

    // Coordinates of all points (one row per point, in ID order) or of the given points
    Matrix getPointCoords() const;
    Matrix getPointCoords(const std::vector<int> &IDs) const;

    // min_x, max_x, min_y, max_y
    std::tuple<double, double, double, double> bounding_box_2d() const;

    // Points on the given curves; points and curves on the given surfaces (including holes);
    // points, curves and surfaces of the given volumes
    std::vector<int> pointsOnCurves(const std::vector<int> &IDs) const;
    std::tuple<std::vector<int>, std::vector<int>> stuffOnSurfaces(const std::vector<int> &IDs) const;
    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> stuffOnVolumes(const std::vector<int> &IDs) const;

private:
    void addCurve(CurveType type, const std::vector<int> &points, int ID, int marker, int el_on_curve,
                  const std::string &el_distrib_type, double el_distrib_val);
    void addSurface(const std::string &name, const std::vector<int> &outer_loop,
                    const std::vector<std::vector<int>> &holes, int ID, int marker, bool is_structured);
    void addVolume(const std::vector<int> &outer_surfaces, const std::vector<std::vector<int>> &holes, int ID,
                   int marker, bool is_structured);

    template <typename T> static int smallestFreeKey(const std::map<int, T> &m);

    std::map<int, Point> m_points;
    std::map<int, Curve> m_curves;
    std::map<int, Surface> m_surfaces;
    std::map<int, Volume> m_volumes;

    bool m_is3D = false;

    bool m_pointIDspecified = false, m_curveIDspecified = false;
    bool m_surfaceIDspecified = false, m_volumeIDspecified = false;
    int m_nextPointID = 0, m_nextCurveID = 0, m_nextSurfaceID = 0, m_nextVolumeID = 0;
};

} // namespace calfem
