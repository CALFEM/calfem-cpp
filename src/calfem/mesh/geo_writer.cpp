//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Gmsh .geo writer, ported from calfem-python (mesh.py, GmshMeshGenerator._writeGeoFile).

#include <calfem/mesh.h>

#include <Eigen/Dense>

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace calfem {

namespace {

const std::set<int> elementsWithQuadFaces = {3, 5, 10, 12, 16, 17, 92, 93};
const std::set<int> secondOrderIncomplete = {9, 11, 13, 14, 16, 17, 18, 19};

std::string number(double v)
{
    char buf[64];
    auto [end, ec] = std::to_chars(buf, buf + sizeof(buf), v);
    return std::string(buf, end);
}

// Comma-separated list; offset is added to positive and subtracted from negative values
std::string formatList(const std::vector<int> &values, int offset = 0)
{
    std::string s;
    for (size_t i = 0; i < values.size(); i++)
    {
        if (i > 0)
            s += ", ";
        int v = values[i];
        s += std::to_string(v < 0 ? v - offset : v + offset);
    }
    return s;
}

const char *curveName(Geometry::CurveType type)
{
    switch (type)
    {
    case Geometry::CurveType::Spline:
        return "Spline";
    case Geometry::CurveType::BSpline:
        return "BSpline";
    case Geometry::CurveType::Circle:
        return "Circle";
    case Geometry::CurveType::Ellipse:
        return "Ellipse";
    }
    return "Spline";
}

// Markers in order of first appearance (calfem-python uses insertion-ordered dicts)
class MarkerGroups {
public:
    void add(int marker, int id)
    {
        if (m_ids.find(marker) == m_ids.end())
            m_order.push_back(marker);
        m_ids[marker].push_back(id);
    }

    template <typename F> void forEach(F f) const
    {
        for (int marker : m_order)
            f(marker, m_ids.at(marker));
    }

private:
    std::vector<int> m_order;
    std::map<int, std::vector<int>> m_ids;
};

class GeoWriter {
public:
    GeoWriter(const Geometry &g, int elType) : m_g(g), m_elType(elType)
    {}

    std::string write()
    {
        MarkerGroups pointMarkers, curveMarkers, surfaceMarkers, volumeMarkers;

        for (const auto &[ID, p] : m_g.points())
        {
            m_out << "Point(" << ID + 1 << ") = {" << number(p.coord(0)) << ", " << number(p.coord(1)) << ", "
                  << number(p.coord(2)) << ", " << number(p.el_size) << "};\n";
            pointMarkers.add(p.marker, ID);
        }

        for (const auto &[ID, c] : m_g.curves())
        {
            m_out << curveName(c.type) << "(" << ID + 1 << ") = {" << formatList(c.points, 1) << "};\n";

            if (c.el_on_curve > 0)
            {
                std::string distribution;
                if (!c.el_distrib_type.empty())
                {
                    char buf[64];
                    std::snprintf(buf, sizeof(buf), "%f", c.el_distrib_val);
                    distribution = "Using " + c.el_distrib_type + " " + buf;
                }
                m_out << "Transfinite Line{" << ID + 1 << "} = " << c.el_on_curve + 1 << " " << distribution << ";\n";
            }

            curveMarkers.add(c.marker, ID);
        }

        for (const auto &[ID, s] : m_g.surfaces())
        {
            writeLineLoop(s.outer_loop, ID + 1);

            std::vector<int> holeIDs;
            for (size_t i = 0; i < s.holes.size(); i++)
            {
                int holeID = 10000 * (ID + 1) + 10 * static_cast<int>(i) + 5;
                writeLineLoop(s.holes[i], holeID);
                holeIDs.push_back(holeID);
            }

            std::string holeString = holeIDs.empty() ? "" : ", " + formatList(holeIDs);
            m_out << s.name << "(" << ID + 1 << ") = {" << ID + 1 << holeString << "};\n";

            if (s.is_structured)
            {
                std::set<int> corners;
                for (int c : s.outer_loop)
                {
                    const auto &pts = m_g.curves().at(c).points;
                    corners.insert(pts.front());
                    corners.insert(pts.back());
                }
                m_out << "Transfinite Surface{" << ID + 1 << "} = {"
                      << formatList(std::vector<int>(corners.begin(), corners.end()), 1) << "};\n";
            }

            surfaceMarkers.add(s.marker, ID);
        }

        for (const auto &[ID, v] : m_g.volumes())
        {
            m_out << "Surface Loop(" << ID + 1 << ") = {" << formatList(v.outer_surfaces, 1) << "};\n";

            std::vector<int> holeIDs;
            for (size_t i = 0; i < v.holes.size(); i++)
            {
                int holeID = 10000 * (ID + 1) + 10 * static_cast<int>(i) + 7;
                m_out << "Surface Loop(" << holeID << ") = {" << formatList(v.holes[i], 1) << "};\n";
                holeIDs.push_back(holeID);
            }

            std::string holeString = holeIDs.empty() ? "" : " , " + formatList(holeIDs);
            m_out << "Volume(" << ID + 1 << ") = {" << ID + 1 << holeString << "};\n";

            if (v.is_structured)
                m_out << "Transfinite Volume{" << ID + 1 << "} = {};\n";

            volumeMarkers.add(v.marker, ID);
        }

        if (elementsWithQuadFaces.count(m_elType))
            m_out << "Mesh.RecombineAll = 1;\n";

        pointMarkers.forEach([this](int marker, const std::vector<int> &ids) {
            if (marker != 0)
                m_out << "Physical Point(" << marker << ") = {" << formatList(ids, 1) << "};\n";
        });
        curveMarkers.forEach([this](int marker, const std::vector<int> &ids) {
            m_out << "Physical Line(" << marker << ") = {" << formatList(ids, 1) << "};\n";
        });
        surfaceMarkers.forEach([this](int marker, const std::vector<int> &ids) {
            m_out << "Physical Surface(" << marker << ") = {" << formatList(ids, 1) << "};\n";
        });
        volumeMarkers.forEach([this](int marker, const std::vector<int> &ids) {
            m_out << "Physical Volume(" << marker << ") = {" << formatList(ids, 1) << "};\n";
        });

        if (secondOrderIncomplete.count(m_elType))
            m_out << "Mesh.SecondOrderIncomplete=1;\n";

        return m_out.str();
    }

private:
    const Geometry::Curve &curve(int id) const
    {
        auto it = m_g.curves().find(id);
        if (it == m_g.curves().end())
            throw std::invalid_argument("Mesh: curve " + std::to_string(id) + " used in a loop does not exist.");
        return it->second;
    }

    Eigen::Vector3d coord(int pointID) const
    {
        return m_g.points().at(pointID).coord.transpose();
    }

    void writeLineLoop(const std::vector<int> &curves, int loopID)
    {
        std::vector<std::vector<int>> endPoints;
        for (int c : curves)
            endPoints.push_back({curve(c).points.front(), curve(c).points.back()});

        std::vector<int> lineIndices;
        for (int c : curves)
            lineIndices.push_back(c < 0 ? c - 1 : c + 1);

        auto contains = [](const std::vector<int> &v, int x) { return std::find(v.begin(), v.end(), x) != v.end(); };

        auto n = curves.size();
        bool isFirstLine = true;

        for (size_t k = 0; k < n; k++)
        {
            if (isFirstLine && (n > 1))
            {
                isFirstLine = false;
                if (contains(endPoints[1], endPoints[0][1]))
                    ;
                else if (contains(endPoints[1], endPoints[0][0]))
                {
                    std::reverse(endPoints[0].begin(), endPoints[0].end());
                    lineIndices[0] *= -1;
                }
                else
                    throw std::invalid_argument("Mesh: the first curve of line loop " + std::to_string(loopID) +
                                                " does not link up to the subsequent curve.");
            }
            else if (endPoints[k][0] == endPoints[k - 1][1])
                ;
            else if (endPoints[k][1] == endPoints[k - 1][1])
            {
                std::reverse(endPoints[k].begin(), endPoints[k].end());
                lineIndices[k] *= -1;
            }
            else
                throw std::invalid_argument("Mesh: curve " + std::to_string(k) + " of line loop " +
                                            std::to_string(loopID) + " does not link up with the preceding curve.");

            if ((k == n - 1) && (endPoints[k][1] != endPoints[0][0]))
                throw std::invalid_argument("Mesh: the last curve of line loop " + std::to_string(loopID) +
                                            " does not join up with the first curve.");
        }

        if (!m_g.is3D())
            lineIndices = makeCounterClockwise(lineIndices);

        m_out << "Line Loop(" << loopID << ") = {" << formatList(lineIndices) << "};\n";
    }

    // Negates all indices if the loop is clockwise (as calfem-python's _makeCounterClockwise)
    std::vector<int> makeCounterClockwise(std::vector<int> lineIndices) const
    {
        double summa = 0.0;

        auto addSegments = [&summa](const std::vector<Eigen::Vector3d> &pts) {
            for (size_t i = 0; i + 1 < pts.size(); i++)
                summa += (pts[i + 1](0) - pts[i](0)) * (pts[i + 1](1) + pts[i](1));
        };

        for (int index : lineIndices)
        {
            int sign = index < 0 ? -1 : 1;
            int realIndex = sign * index - 1;
            const auto &c = curve(realIndex);

            std::vector<Eigen::Vector3d> pts;
            if ((c.type == Geometry::CurveType::Spline) || (c.type == Geometry::CurveType::BSpline))
            {
                for (int id : c.points)
                    pts.push_back(coord(id));
            }
            else
            {
                // Arc approximated by start - "middle" - end, as calfem-python
                Eigen::Vector3d a = coord(c.points[0]), b = coord(c.points[1]);
                Eigen::Vector3d e = coord(c.type == Geometry::CurveType::Circle ? c.points[2] : c.points[3]);
                double r = c.type == Geometry::CurveType::Circle ? (a - b).norm()
                                                                  : ((a - b).norm() + (e - b).norm()) / 2.0;
                Eigen::Vector3d s = a + 2.0 * b + e;
                Eigen::Vector3d d = b + r * s / s.norm();
                pts = {a, d, e};
            }

            if (sign == -1)
                std::reverse(pts.begin(), pts.end());
            addSegments(pts);
        }

        if (summa > 0.0)
            for (auto &x : lineIndices)
                x = -x;

        return lineIndices;
    }

    const Geometry &m_g;
    int m_elType;
    std::ostringstream m_out;
};

} // namespace

std::string geo_string(const Geometry &geometry, int el_type)
{
    return GeoWriter(geometry, el_type).write();
}

} // namespace calfem
