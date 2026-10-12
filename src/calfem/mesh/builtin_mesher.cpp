//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Built-in 2D mesh generator: structured (transfinite) surfaces and constrained Delaunay
// triangulation (CDT library) of general surfaces with graded interior points.

#include "curves.h"
#include "mesh_internal.h"

#include <CDT.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calfem::detail {

namespace {

using Vec2 = Eigen::Vector2d;

// ---- Node registry ----

class NodeRegistry {
public:
    int add(const Vec3 &p)
    {
        m_nodes.push_back(p);
        return static_cast<int>(m_nodes.size()) - 1;
    }

    const Vec3 &operator[](int i) const
    {
        return m_nodes[i];
    }

    Vec3 &operator[](int i)
    {
        return m_nodes[i];
    }

    Vec2 xy(int i) const
    {
        return {m_nodes[i](0), m_nodes[i](1)};
    }

    std::vector<Vec3> &nodes()
    {
        return m_nodes;
    }

private:
    std::vector<Vec3> m_nodes;
};

double cross2(const Vec2 &a, const Vec2 &b)
{
    return a(0) * b(1) - a(1) * b(0);
}

double signedArea(const NodeRegistry &reg, const std::vector<int> &polygon)
{
    double a = 0.0;
    for (size_t i = 0; i < polygon.size(); i++)
        a += cross2(reg.xy(polygon[i]), reg.xy(polygon[(i + 1) % polygon.size()]));
    return 0.5 * a;
}

double triangleArea(const NodeRegistry &reg, int a, int b, int c)
{
    return 0.5 * cross2(reg.xy(b) - reg.xy(a), reg.xy(c) - reg.xy(a));
}

// ---- Curve loops ----

struct OrientedCurve {
    int id;
    bool reversed;
};

// Orders the curves of a loop head to tail (as calfem-python's _writeLineLoop)
std::vector<OrientedCurve> orientLoop(const Geometry &g, const std::vector<int> &loop, int surfaceID)
{
    auto ends = [&](int id) {
        auto it = g.curves().find(id);
        if (it == g.curves().end())
            throw std::invalid_argument("Mesh: surface " + std::to_string(surfaceID) + " uses curve " +
                                        std::to_string(id) + ", which does not exist.");
        return std::make_pair(it->second.points.front(), it->second.points.back());
    };

    std::vector<OrientedCurve> result;
    std::vector<std::pair<int, int>> e;
    for (int id : loop)
        e.push_back(ends(id));

    auto n = loop.size();
    std::vector<bool> reversed(n, false);

    for (size_t k = 0; k < n; k++)
    {
        if ((k == 0) && (n > 1))
        {
            if ((e[0].second == e[1].first) || (e[0].second == e[1].second))
                ;
            else if ((e[0].first == e[1].first) || (e[0].first == e[1].second))
            {
                std::swap(e[0].first, e[0].second);
                reversed[0] = true;
            }
            else
                throw std::invalid_argument("Mesh: the first curve of a loop of surface " + std::to_string(surfaceID) +
                                            " does not connect to the next curve.");
        }
        else if (k > 0)
        {
            if (e[k].first == e[k - 1].second)
                ;
            else if (e[k].second == e[k - 1].second)
            {
                std::swap(e[k].first, e[k].second);
                reversed[k] = true;
            }
            else
                throw std::invalid_argument("Mesh: curve " + std::to_string(loop[k]) + " in a loop of surface " +
                                            std::to_string(surfaceID) + " does not connect to the previous curve.");
        }

        if ((k == n - 1) && (e[k].second != e[0].first))
            throw std::invalid_argument("Mesh: the loop of surface " + std::to_string(surfaceID) + " is not closed.");
    }

    for (size_t k = 0; k < n; k++)
        result.push_back({loop[k], reversed[k]});
    return result;
}

// ---- Size field for unstructured meshing ----

class SizeField {
public:
    SizeField(const std::vector<Vec2> &points, const std::vector<double> &sizes) : m_points(points), m_sizes(sizes)
    {
        m_hmax = *std::max_element(sizes.begin(), sizes.end());
        m_hmin = *std::min_element(sizes.begin(), sizes.end());
    }

    // Grows from the boundary sizes with the gradation rate, limited to the largest boundary size
    double operator()(const Vec2 &x) const
    {
        double h = m_hmax;
        for (size_t i = 0; i < m_points.size(); i++)
            h = std::min(h, m_sizes[i] + gradation * (x - m_points[i]).norm());
        return h;
    }

    double hmin() const
    {
        return m_hmin;
    }

    static constexpr double gradation = 0.3;

private:
    std::vector<Vec2> m_points;
    std::vector<double> m_sizes;
    double m_hmax = 0.0, m_hmin = 0.0;
};

bool pointInPolygon(const Vec2 &p, const std::vector<Vec2> &poly)
{
    bool inside = false;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
    {
        const Vec2 &a = poly[i], &b = poly[j];
        if (((a(1) > p(1)) != (b(1) > p(1))) && (p(0) < (b(0) - a(0)) * (p(1) - a(1)) / (b(1) - a(1)) + a(0)))
            inside = !inside;
    }
    return inside;
}

double distanceToSegment(const Vec2 &p, const Vec2 &a, const Vec2 &b)
{
    Vec2 d = b - a;
    double t = std::clamp((p - a).dot(d) / d.squaredNorm(), 0.0, 1.0);
    return (p - (a + t * d)).norm();
}

// Deterministic uniform numbers (std::mt19937 is specified exactly; distributions are not)
class Random {
public:
    double uniform()
    {
        return static_cast<double>(m_gen()) / 4294967296.0;
    }

    size_t index(size_t n)
    {
        return static_cast<size_t>(uniform() * static_cast<double>(n)) % n;
    }

private:
    std::mt19937 m_gen{20261012u};
};

// Interior points with spacing following the size field (Poisson disk sampling, seeded with
// the boundary nodes)
std::vector<Vec2> interiorPoints(const std::vector<std::vector<Vec2>> &loops, const SizeField &h)
{
    const double radiusFactor = 0.85;

    Vec2 lo = loops[0][0], hi = loops[0][0];
    for (const auto &loop : loops)
        for (const auto &p : loop)
        {
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }

    double cell = std::max(h.hmin() * radiusFactor, (hi - lo).norm() / 2000.0);
    auto nx = static_cast<int>((hi(0) - lo(0)) / cell) + 1;
    auto ny = static_cast<int>((hi(1) - lo(1)) / cell) + 1;
    std::vector<std::vector<int>> grid(static_cast<size_t>(nx) * ny);

    std::vector<Vec2> pts;
    std::vector<double> radius;

    auto cellOf = [&](const Vec2 &p) {
        int i = std::clamp(static_cast<int>((p(0) - lo(0)) / cell), 0, nx - 1);
        int j = std::clamp(static_cast<int>((p(1) - lo(1)) / cell), 0, ny - 1);
        return std::make_pair(i, j);
    };

    auto insert = [&](const Vec2 &p, double r) {
        pts.push_back(p);
        radius.push_back(r);
        auto [i, j] = cellOf(p);
        grid[static_cast<size_t>(j) * nx + i].push_back(static_cast<int>(pts.size()) - 1);
    };

    auto inside = [&](const Vec2 &p) {
        if (!pointInPolygon(p, loops[0]))
            return false;
        for (size_t k = 1; k < loops.size(); k++)
            if (pointInPolygon(p, loops[k]))
                return false;
        return true;
    };

    // Far enough from all existing points (using the smaller of the two radii)
    auto acceptable = [&](const Vec2 &p, double r) {
        int reach = static_cast<int>(std::ceil(r / cell));
        auto [ci, cj] = cellOf(p);
        for (int j = std::max(0, cj - reach); j <= std::min(ny - 1, cj + reach); j++)
            for (int i = std::max(0, ci - reach); i <= std::min(nx - 1, ci + reach); i++)
                for (int q : grid[static_cast<size_t>(j) * nx + i])
                    if ((pts[q] - p).norm() < std::min(r, radius[q]))
                        return false;
        return true;
    };

    // Seed with the boundary
    for (const auto &loop : loops)
        for (const auto &p : loop)
            insert(p, radiusFactor * h(p));

    auto nBoundary = pts.size();

    // Distance from the boundary segments, to avoid flat triangles at the boundary
    auto farFromBoundary = [&](const Vec2 &p, double r) {
        for (const auto &loop : loops)
            for (size_t i = 0; i < loop.size(); i++)
            {
                const Vec2 &a = loop[i], &b = loop[(i + 1) % loop.size()];
                if ((p - a).norm() > (b - a).norm() + r)
                    continue;
                if (distanceToSegment(p, a, b) < 0.5 * r)
                    return false;
            }
        return true;
    };

    Random rnd;
    std::vector<int> active;
    for (size_t i = 0; i < pts.size(); i++)
        active.push_back(static_cast<int>(i));

    const int attempts = 30;
    const double twoPi = 6.28318530717958647692;

    while (!active.empty())
    {
        size_t k = rnd.index(active.size());
        int a = active[k];
        bool found = false;

        for (int t = 0; t < attempts; t++)
        {
            double ra = radius[a];
            double dist = ra * (1.0 + rnd.uniform());
            double angle = twoPi * rnd.uniform();
            Vec2 c = pts[a] + dist * Vec2(std::cos(angle), std::sin(angle));

            if (!inside(c))
                continue;

            double rc = radiusFactor * h(c);
            if (!acceptable(c, rc) || !farFromBoundary(c, rc))
                continue;

            insert(c, rc);
            active.push_back(static_cast<int>(pts.size()) - 1);
            found = true;
            break;
        }

        if (!found)
        {
            active[k] = active.back();
            active.pop_back();
        }
    }

    return std::vector<Vec2>(pts.begin() + static_cast<std::ptrdiff_t>(nBoundary), pts.end());
}

// Constrained Delaunay triangulation of the nodes (global indices) with the loop edges as
// constraints; returns triangles (counter-clockwise) in global node indices
std::vector<std::array<int, 3>> triangulate(const NodeRegistry &reg, const std::vector<int> &nodes,
                                            const std::vector<std::pair<int, int>> &edges)
{
    std::map<int, CDT::VertInd> local;
    std::vector<CDT::V2d<double>> vertices;
    for (int n : nodes)
    {
        local[n] = static_cast<CDT::VertInd>(vertices.size());
        vertices.emplace_back(reg[n](0), reg[n](1));
    }

    std::vector<CDT::Edge> cdtEdges;
    for (const auto &[a, b] : edges)
        cdtEdges.emplace_back(local.at(a), local.at(b));

    CDT::Triangulation<double> cdt;
    cdt.insertVertices(vertices);
    cdt.insertEdges(cdtEdges);
    cdt.eraseOuterTrianglesAndHoles();

    std::vector<std::array<int, 3>> triangles;
    for (const auto &t : cdt.triangles)
        triangles.push_back({nodes[t.vertices[0]], nodes[t.vertices[1]], nodes[t.vertices[2]]});
    return triangles;
}

// Laplacian smoothing of the free nodes; a move is kept only if all adjacent elements stay valid
template <typename Elements, typename Valid>
void smooth(NodeRegistry &reg, const Elements &elements, const std::set<int> &fixed, int iterations,
            const Valid &valid)
{
    std::map<int, std::set<int>> neighbours;
    std::map<int, std::vector<size_t>> adjacent;

    for (size_t e = 0; e < elements.size(); e++)
    {
        const auto &el = elements[e];
        auto n = el.size();
        for (size_t i = 0; i < n; i++)
        {
            adjacent[el[i]].push_back(e);
            neighbours[el[i]].insert(el[(i + 1) % n]);
            neighbours[el[i]].insert(el[(i + n - 1) % n]);
        }
    }

    for (int it = 0; it < iterations; it++)
        for (auto &[node, nb] : neighbours)
        {
            if (fixed.count(node))
                continue;

            Vec3 old = reg[node];
            Vec3 avg = Vec3::Zero();
            for (int q : nb)
                avg += reg[q];
            reg[node] = avg / static_cast<double>(nb.size());

            for (size_t e : adjacent[node])
                if (!valid(elements[e]))
                {
                    reg[node] = old;
                    break;
                }
        }
}

} // namespace

RawMesh meshBuiltIn(const Geometry &g, int elType, double factor)
{
    if (g.is3D() || !g.volumes().empty())
        throw std::invalid_argument("Mesh: the built-in mesher supports 2D geometries only; use the Gmsh backend "
                                    "for 3D (build calfem with CALFEM_WITH_GMSH).");

    if ((elType != 2) && (elType != 3))
        throw std::invalid_argument("Mesh: the built-in mesher supports el_type 2 (triangles) and 3 "
                                    "(quadrilaterals); other element types require the Gmsh backend.");

    if (g.surfaces().empty())
        throw std::invalid_argument("Mesh: the geometry has no surfaces.");

    RawMesh raw;
    raw.dim = 2;
    NodeRegistry reg;

    // Curves that must be divided into an even number of elements (quads on unstructured surfaces)
    std::set<int> evenCurves;
    if (elType == 3)
        for (const auto &[id, s] : g.surfaces())
            if (!s.is_structured)
            {
                evenCurves.insert(s.outer_loop.begin(), s.outer_loop.end());
                for (const auto &hole : s.holes)
                    evenCurves.insert(hole.begin(), hole.end());
            }

    // ---- Points (curve end points) ----

    std::map<int, int> pointNode;
    std::set<int> curvePoints;
    for (const auto &[id, c] : g.curves())
    {
        curvePoints.insert(c.points.front());
        curvePoints.insert(c.points.back());
    }

    for (int pid : curvePoints)
    {
        auto it = g.points().find(pid);
        if (it == g.points().end())
            throw std::invalid_argument("Mesh: point " + std::to_string(pid) + " used by a curve does not exist.");
        pointNode[pid] = reg.add(it->second.coord.transpose());
    }

    // ---- Curves ----

    std::map<int, std::vector<int>> curveNodes;

    for (const auto &[id, c] : g.curves())
    {
        CurveShape shape(g, id);

        double h0 = g.points().at(c.points.front()).el_size * factor;
        double h1 = g.points().at(c.points.back()).el_size * factor;
        bool even = evenCurves.count(id) > 0;

        auto xi = curveDivision(c, shape.length(), h0, h1, even);
        auto n = static_cast<int>(xi.size()) - 1;

        if (even && (n % 2 == 1))
            throw std::invalid_argument("Mesh: curve " + std::to_string(id) +
                                        " needs an even el_on_curve for quadrilaterals (el_type 3) on an "
                                        "unstructured surface.");

        std::vector<int> nodes;
        nodes.push_back(pointNode.at(c.points.front()));
        for (int k = 1; k < n; k++)
            nodes.push_back(reg.add(shape.atLength(xi[k])));
        nodes.push_back(pointNode.at(c.points.back()));
        curveNodes[id] = nodes;

        // Boundary nodes and line elements (every curve, also marker 0, as calfem-python)
        for (int node : nodes)
            raw.boundaryNodes[c.marker].insert(node + 1);
        for (size_t k = 0; k + 1 < nodes.size(); k++)
            raw.boundaryElements[c.marker].push_back({1, {nodes[k] + 1, nodes[k + 1] + 1}});
    }

    // Marked points
    for (const auto &[pid, node] : pointNode)
    {
        int marker = g.points().at(pid).marker;
        if (marker != 0)
        {
            raw.boundaryNodes[marker].insert(node + 1);
            raw.boundaryElements[marker].push_back({15, {node + 1}});
        }
    }

    // Nodes along an oriented loop (closed, last node not repeated)
    auto loopNodes = [&](const std::vector<OrientedCurve> &loop) {
        std::vector<int> nodes;
        for (const auto &oc : loop)
        {
            std::vector<int> cn = curveNodes.at(oc.id);
            if (oc.reversed)
                std::reverse(cn.begin(), cn.end());
            nodes.insert(nodes.end(), cn.begin(), cn.end() - 1);
        }
        return nodes;
    };

    auto addElement = [&](std::vector<int> nodes, int marker) {
        for (auto &n : nodes)
            n += 1;
        raw.elements.push_back(nodes);
        raw.markers.push_back(marker);
    };

    // ---- Surfaces ----

    for (const auto &[sid, s] : g.surfaces())
    {
        auto loop = orientLoop(g, s.outer_loop, sid);
        auto outer = loopNodes(loop);

        // Counter-clockwise outer boundary
        if (signedArea(reg, outer) < 0.0)
        {
            std::reverse(loop.begin(), loop.end());
            for (auto &oc : loop)
                oc.reversed = !oc.reversed;
            outer = loopNodes(loop);
        }

        if (s.is_structured)
        {
            // ---- Transfinite (mapped) mesh ----

            std::vector<std::vector<int>> side(4);
            for (int k = 0; k < 4; k++)
            {
                side[k] = curveNodes.at(loop[k].id);
                if (loop[k].reversed)
                    std::reverse(side[k].begin(), side[k].end());
            }

            auto n0 = static_cast<int>(side[0].size()) - 1;
            auto n1 = static_cast<int>(side[1].size()) - 1;

            if ((static_cast<int>(side[2].size()) - 1 != n0) || (static_cast<int>(side[3].size()) - 1 != n1))
                throw std::invalid_argument("Mesh: opposite curves of structured surface " + std::to_string(sid) +
                                            " have different numbers of elements.");

            auto B = [&](int i) { return side[0][i]; };
            auto R = [&](int j) { return side[1][j]; };
            auto T = [&](int i) { return side[2][n0 - i]; };
            auto L = [&](int j) { return side[3][n1 - j]; };

            // Arc length fractions along the sides
            auto fractions = [&](auto node, int n) {
                std::vector<double> f(n + 1, 0.0);
                for (int i = 1; i <= n; i++)
                    f[i] = f[i - 1] + (reg[node(i)] - reg[node(i - 1)]).norm();
                for (auto &v : f)
                    v /= f[n];
                return f;
            };

            auto ub = fractions(B, n0), ut = fractions(T, n0);
            auto vl = fractions(L, n1), vr = fractions(R, n1);

            Vec3 P00 = reg[B(0)], P10 = reg[B(n0)], P01 = reg[T(0)], P11 = reg[T(n0)];

            std::vector<std::vector<int>> grid(n0 + 1, std::vector<int>(n1 + 1));
            for (int i = 0; i <= n0; i++)
                for (int j = 0; j <= n1; j++)
                {
                    if (j == 0)
                        grid[i][j] = B(i);
                    else if (j == n1)
                        grid[i][j] = T(i);
                    else if (i == 0)
                        grid[i][j] = L(j);
                    else if (i == n0)
                        grid[i][j] = R(j);
                    else
                    {
                        // Coupled parameters for graded sides
                        double u = 0.5 * (ub[i] + ut[i]), v = 0.5 * (vl[j] + vr[j]);
                        for (int iter = 0; iter < 4; iter++)
                        {
                            u = (1.0 - v) * ub[i] + v * ut[i];
                            v = (1.0 - u) * vl[j] + u * vr[j];
                        }

                        Vec3 p = (1.0 - v) * reg[B(i)] + v * reg[T(i)] + (1.0 - u) * reg[L(j)] + u * reg[R(j)] -
                                 ((1.0 - u) * (1.0 - v) * P00 + u * (1.0 - v) * P10 + (1.0 - u) * v * P01 +
                                  u * v * P11);
                        grid[i][j] = reg.add(p);
                    }
                }

            for (int j = 0; j < n1; j++)
                for (int i = 0; i < n0; i++)
                {
                    int a = grid[i][j], b = grid[i + 1][j], c = grid[i + 1][j + 1], d = grid[i][j + 1];
                    if (elType == 3)
                        addElement({a, b, c, d}, s.marker);
                    else
                    {
                        addElement({a, b, c}, s.marker);
                        addElement({a, c, d}, s.marker);
                    }
                }

            continue;
        }

        // ---- Unstructured mesh ----

        std::vector<std::vector<int>> loops{outer};
        for (const auto &hole : s.holes)
            loops.push_back(loopNodes(orientLoop(g, hole, sid)));

        // For quadrilaterals: triangulate the coarse boundary (every other node) and split each
        // triangle into three quadrilaterals using the remaining boundary nodes as edge midpoints
        std::map<std::pair<int, int>, int> boundaryMidpoint;
        std::vector<std::vector<int>> coarse = loops;

        if (elType == 3)
            for (auto &ring : coarse)
            {
                std::vector<int> c;
                for (size_t i = 0; i < ring.size(); i += 2)
                {
                    int a = ring[i], m = ring[i + 1], b = ring[(i + 2) % ring.size()];
                    boundaryMidpoint[{std::min(a, b), std::max(a, b)}] = m;
                    c.push_back(a);
                }
                ring = c;
            }

        // Size field from the boundary spacing
        std::vector<Vec2> bpoints;
        std::vector<double> bsizes;
        std::vector<std::vector<Vec2>> polygons;
        std::set<int> fixed;

        for (const auto &ring : coarse)
        {
            std::vector<Vec2> poly;
            for (size_t i = 0; i < ring.size(); i++)
            {
                int prev = ring[(i + ring.size() - 1) % ring.size()], next = ring[(i + 1) % ring.size()];
                bpoints.push_back(reg.xy(ring[i]));
                bsizes.push_back(0.5 * ((reg.xy(ring[i]) - reg.xy(prev)).norm() +
                                        (reg.xy(ring[i]) - reg.xy(next)).norm()));
                poly.push_back(reg.xy(ring[i]));
                fixed.insert(ring[i]);
            }
            polygons.push_back(poly);
        }

        SizeField h(bpoints, bsizes);

        std::vector<int> nodes;
        std::vector<std::pair<int, int>> edges;
        for (const auto &ring : coarse)
            for (size_t i = 0; i < ring.size(); i++)
            {
                nodes.push_back(ring[i]);
                edges.emplace_back(ring[i], ring[(i + 1) % ring.size()]);
            }

        for (const auto &p : interiorPoints(polygons, h))
            nodes.push_back(reg.add(Vec3(p(0), p(1), 0.0)));

        auto positive = [&](const std::array<int, 3> &t) { return triangleArea(reg, t[0], t[1], t[2]) > 0.0; };

        // Triangulate, smooth, re-triangulate (restores the Delaunay property), smooth again
        auto triangles = triangulate(reg, nodes, edges);
        smooth(reg, triangles, fixed, 5, positive);
        triangles = triangulate(reg, nodes, edges);
        smooth(reg, triangles, fixed, 3, positive);

        if (elType == 2)
        {
            for (const auto &t : triangles)
                addElement({t[0], t[1], t[2]}, s.marker);
            continue;
        }

        // Split into quadrilaterals
        std::map<std::pair<int, int>, int> midpoint = boundaryMidpoint;
        auto mid = [&](int a, int b) {
            std::pair<int, int> key{std::min(a, b), std::max(a, b)};
            auto it = midpoint.find(key);
            if (it != midpoint.end())
                return it->second;
            int m = reg.add(0.5 * (reg[a] + reg[b]));
            midpoint[key] = m;
            return m;
        };

        std::vector<std::array<int, 4>> quads;
        for (const auto &t : triangles)
        {
            int a = t[0], b = t[1], c = t[2];
            int g0 = reg.add((reg[a] + reg[b] + reg[c]) / 3.0);
            int mab = mid(a, b), mbc = mid(b, c), mca = mid(c, a);
            quads.push_back({a, mab, g0, mca});
            quads.push_back({b, mbc, g0, mab});
            quads.push_back({c, mca, g0, mbc});
        }

        std::set<int> fixedQuad;
        for (const auto &ring : loops)
            fixedQuad.insert(ring.begin(), ring.end());

        auto positiveQuad = [&](const std::array<int, 4> &q) {
            return (triangleArea(reg, q[0], q[1], q[2]) > 0.0) && (triangleArea(reg, q[0], q[2], q[3]) > 0.0) &&
                   (triangleArea(reg, q[0], q[1], q[3]) > 0.0) && (triangleArea(reg, q[1], q[2], q[3]) > 0.0);
        };

        smooth(reg, quads, fixedQuad, 10, positiveQuad);

        for (const auto &q : quads)
            addElement({q[0], q[1], q[2], q[3]}, s.marker);
    }

    raw.nodes = reg.nodes();
    return raw;
}

} // namespace calfem::detail
