//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// calfem-python mesh and field drawing functions (draw_mesh, draw_element_values,
// draw_nodal_values_contourf, ..., eliso2, elval2, elflux2), ported from calfem/vis_mpl.py.
// matplotlib's tricontour/tricontourf/tripcolor/quiver are replaced by contouring on the
// triangulated mesh, with the level selection and colour mapping of matplotlib.

#include <calfem/vis.h>

#include "detail/log.h"
#include "detail/vis_codes.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace calfem::vis {

namespace {

const double inf = std::numeric_limits<double>::infinity();

// ---- Mesh helpers ----

void setTitle(Figure &fig, const std::string &title)
{
    if (!title.empty())
        fig.set_title(title);
}

// Node index (0-based) of dof (1-based, numbered node by node)
int dofNode(int dof, int dofsPerNode)
{
    return (dof - 1) / dofsPerNode;
}

// Faces (0-based node indices) of the elements, as calfem-python's ce2vf for 2D elements
std::vector<std::vector<int>> meshFaces(const IntMatrix &edof, int dofsPerNode, int elType, Eigen::Index nNodes,
                                        const char *caller)
{
    int faceNodes;
    switch (elType)
    {
    case 2:
        faceNodes = 3;
        break;
    case 3:
    case 16:
        faceNodes = 4; // 8-node quadrilaterals: the corner nodes
        break;
    default:
        throw std::invalid_argument(std::string(caller) + ": element type " + std::to_string(elType) +
                                    " is not supported (2 triangles, 3 quadrilaterals, 16 8-node quadrilaterals).");
    }
    if (dofsPerNode < 1)
        throw std::invalid_argument(std::string(caller) + ": dofs_per_node must be positive.");
    if (edof.cols() < faceNodes * dofsPerNode)
        throw std::invalid_argument(std::string(caller) + ": edof has too few columns for the element type.");

    std::vector<std::vector<int>> faces(static_cast<std::size_t>(edof.rows()));
    for (Eigen::Index e = 0; e < edof.rows(); e++)
    {
        auto &face = faces[static_cast<std::size_t>(e)];
        for (int k = 0; k < faceNodes; k++)
        {
            int n = dofNode(edof(e, k * dofsPerNode), dofsPerNode);
            if (n < 0 || n >= nNodes)
                throw std::invalid_argument(std::string(caller) + ": edof refers to a node outside coords.");
            face.push_back(n);
        }
    }
    return faces;
}

// Node topology (1-based) of the elements: edof itself for dofsPerNode = 0, as calfem-python
std::vector<std::vector<int>> nodeTopology(const IntMatrix &edof, int dofsPerNode, Eigen::Index nNodes,
                                           const char *caller)
{
    int step = std::max(1, dofsPerNode);
    if (edof.cols() % step != 0)
        throw std::invalid_argument(std::string(caller) + ": edof columns must be a multiple of dofs_per_node.");

    std::vector<std::vector<int>> topo(static_cast<std::size_t>(edof.rows()));
    for (Eigen::Index e = 0; e < edof.rows(); e++)
        for (Eigen::Index k = 0; k < edof.cols(); k += step)
        {
            int n = dofsPerNode > 0 ? dofNode(edof(e, k), dofsPerNode) : edof(e, k) - 1;
            if (n < 0 || n >= nNodes)
                throw std::invalid_argument(std::string(caller) + ": edof refers to a node outside coords.");
            topo[static_cast<std::size_t>(e)].push_back(n);
        }
    return topo;
}

using Triangle = std::array<int, 3>;

// Triangles of 3-, 4- and 8-node elements (calfem-python's topo_to_tri)
std::vector<Triangle> toTriangles(const std::vector<std::vector<int>> &topo, const char *caller)
{
    std::vector<Triangle> tris;
    for (const auto &t : topo)
    {
        switch (t.size())
        {
        case 3:
            tris.push_back({t[0], t[1], t[2]});
            break;
        case 4:
            tris.push_back({t[0], t[1], t[2]});
            tris.push_back({t[2], t[3], t[0]});
            break;
        case 8:
            tris.push_back({t[0], t[4], t[7]});
            tris.push_back({t[4], t[1], t[5]});
            tris.push_back({t[5], t[2], t[6]});
            tris.push_back({t[6], t[3], t[7]});
            tris.push_back({t[4], t[6], t[7]});
            tris.push_back({t[4], t[5], t[6]});
            break;
        default:
            throw std::invalid_argument(std::string(caller) + ": element topology not supported (3, 4 or 8 nodes).");
        }
    }
    return tris;
}

void checkCoords(const Matrix &coords, const char *caller)
{
    if (coords.cols() != 2)
        throw std::invalid_argument(std::string(caller) + ": only 2D meshes are supported (coords must be N x 2).");
}

void checkNodalValues(const ColVec &values, const Matrix &coords, const char *caller)
{
    if (values.size() != coords.rows())
        throw std::invalid_argument(std::string(caller) + ": there must be one value per node.");
}

std::vector<Point> facePoints(const Matrix &coords, const std::vector<int> &face)
{
    std::vector<Point> points;
    for (int n : face)
        points.push_back({coords(n, 0), coords(n, 1)});
    return points;
}

// Mesh edges as closed lines (calfem-python's draw_mesh without filling)
void addMeshLines(Figure &fig, const Matrix &coords, const std::vector<std::vector<int>> &faces, const Color &color)
{
    for (const auto &face : faces)
    {
        Line line;
        line.points = facePoints(coords, face);
        line.color = color;
        line.closed = true;
        fig.add(std::move(line));
    }
}

// Nodal displacements (N x 2) from the global displacement vector or N x 2 nodal values
// (calfem-python's _nodal_displacements)
Matrix nodalDisplacements(const Matrix &a, const Matrix &coords, int dofsPerNode, const char *caller)
{
    Eigen::Index n = coords.rows();
    Matrix u(n, 2);
    if (a.rows() == n && a.cols() == 2)
        return a;
    if ((a.rows() == 1 || a.cols() == 1) && a.size() == n * dofsPerNode && dofsPerNode >= 2)
    {
        for (Eigen::Index i = 0; i < n; i++)
            for (int d = 0; d < 2; d++)
                u(i, d) = a(i * dofsPerNode + d);
        return u;
    }
    if ((a.rows() == 1 || a.cols() == 1) && a.size() == 2 * n)
    {
        for (Eigen::Index i = 0; i < n; i++)
            for (int d = 0; d < 2; d++)
                u(i, d) = a(2 * i + d);
        return u;
    }
    throw std::invalid_argument(std::string(caller) +
                                ": displacements must be the global displacement vector or N x 2 nodal values.");
}

// ---- Contour levels (matplotlib's MaxNLocator and ContourSet._autolev) ----

// Python's float divmod
std::pair<double, double> pyDivmod(double vx, double wx)
{
    double mod = std::fmod(vx, wx);
    double div = (vx - mod) / wx;
    if (mod != 0.0)
    {
        if ((wx < 0.0) != (mod < 0.0))
        {
            mod += wx;
            div -= 1.0;
        }
    }
    else
        mod = std::copysign(0.0, wx);

    double floordiv;
    if (div != 0.0)
    {
        floordiv = std::floor(div);
        if (div - floordiv > 0.5)
            floordiv += 1.0;
    }
    else
        floordiv = std::copysign(0.0, vx / wx);
    return {floordiv, mod};
}

double pyFloorDiv(double vx, double wx)
{
    return pyDivmod(vx, wx).first;
}

// MaxNLocator(nbins, min_n_ticks = 1).tick_values(vmin, vmax) with the default steps
std::vector<double> maxNLocatorTicks(double vmin, double vmax, int nbins)
{
    // mtransforms.nonsingular(vmin, vmax, expander=1e-13, tiny=1e-14)
    const double expander = 1e-13, tiny = 1e-14;
    if (!std::isfinite(vmin) || !std::isfinite(vmax))
    {
        vmin = -expander;
        vmax = expander;
    }
    else
    {
        if (vmax < vmin)
            std::swap(vmin, vmax);
        double maxabs = std::max(std::abs(vmin), std::abs(vmax));
        if (maxabs < (1e6 / tiny) * DBL_MIN)
        {
            vmin = -expander;
            vmax = expander;
        }
        else if (vmax - vmin <= maxabs * tiny)
        {
            if (vmax == 0.0 && vmin == 0.0)
            {
                vmin = -expander;
                vmax = expander;
            }
            else
            {
                vmin -= expander * std::abs(vmin);
                vmax += expander * std::abs(vmax);
            }
        }
    }

    // scale_range(vmin, vmax, nbins)
    double dv = std::abs(vmax - vmin);
    double meanv = 0.5 * (vmax + vmin);
    double offset = 0.0;
    if (!(std::abs(meanv) / dv < 100.0))
        offset = std::copysign(std::pow(10.0, pyFloorDiv(std::log10(std::abs(meanv)), 1.0)), meanv);
    double scale = std::pow(10.0, pyFloorDiv(std::log10(dv / nbins), 1.0));

    // MaxNLocator's default steps, extended as _staircase: 0.1 * steps[:-1], steps, 10 * steps[1]
    double lo = vmin - offset, hi = vmax - offset;
    const double defaultSteps[] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    std::vector<double> steps;
    for (int i = 0; i < 9; i++)
        steps.push_back(0.1 * defaultSteps[i] * scale);
    for (double s : defaultSteps)
        steps.push_back(s * scale);
    steps.push_back(10.0 * defaultSteps[1] * scale);

    double rawStep = (hi - lo) / nbins;
    std::size_t istep = steps.size() - 1;
    for (std::size_t i = 0; i < steps.size(); i++)
        if (steps[i] >= rawStep)
        {
            istep = i;
            break;
        }

    std::vector<double> ticks;
    for (std::size_t s = istep + 1; s-- > 0;)
    {
        double step = steps[s];

        // _Edge_integer(step, offset)
        auto closeto = [&](double ms, double edge) {
            double tol = 1e-10;
            if (std::abs(offset) > 0.0)
            {
                double digits = std::log10(std::abs(offset) / step);
                tol = std::min(0.4999, std::max(1e-10, std::pow(10.0, digits - 12.0)));
            }
            return std::abs(ms - edge) < tol;
        };
        auto le = [&](double x) {
            auto [d, m] = pyDivmod(x, step);
            return closeto(m / step, 1.0) ? d + 1.0 : d;
        };
        auto ge = [&](double x) {
            auto [d, m] = pyDivmod(x, step);
            return closeto(m / step, 0.0) ? d : d + 1.0;
        };

        double bestVmin = pyFloorDiv(lo, step) * step;
        double low = le(lo - bestVmin), high = ge(hi - bestVmin);

        ticks.clear();
        int nticks = 0;
        for (double k = low; k <= high; k += 1.0)
        {
            double t = k * step + bestVmin;
            ticks.push_back(t);
            if (t <= hi && t >= lo)
                nticks++;
        }
        if (nticks >= 1)
            break;
    }

    for (double &t : ticks)
        t += offset;
    return ticks;
}

// Levels for n contour levels over [zmin, zmax] (ContourSet._autolev(n), extend = "neither")
std::vector<double> autoLevels(double zmin, double zmax, int n)
{
    if (n < 1)
        throw std::invalid_argument("The number of contour levels must be positive.");

    std::vector<double> lev = maxNLocatorTicks(zmin, zmax, n + 1);

    std::size_t i0 = 0, i1 = lev.size();
    for (std::size_t i = 0; i < lev.size(); i++)
        if (lev[i] < zmin)
            i0 = i;
    for (std::size_t i = 0; i < lev.size(); i++)
        if (lev[i] > zmax)
        {
            i1 = i + 1;
            break;
        }
    if (i1 - i0 < 3)
    {
        i0 = 0;
        i1 = lev.size();
    }
    return std::vector<double>(lev.begin() + std::ptrdiff_t(i0), lev.begin() + std::ptrdiff_t(i1));
}

void checkLevels(const std::vector<double> &levels, bool filled, const char *caller)
{
    if (levels.size() < (filled ? 2u : 1u))
        throw std::invalid_argument(std::string(caller) +
                                    (filled ? ": filled contours require at least 2 levels." : ": no contour levels."));
    for (std::size_t i = 1; i < levels.size(); i++)
        if (!(levels[i] > levels[i - 1]))
            throw std::invalid_argument(std::string(caller) + ": contour levels must be increasing.");
}

std::pair<double, double> valueRange(const ColVec &values)
{
    if (values.size() == 0)
        return {0.0, 0.0};
    return {values.minCoeff(), values.maxCoeff()};
}

// ---- Filled contours on triangles ----

struct ValuePoint {
    double x;
    double y;
    double v;
};

// Part of the (convex) polygon where v >= level (keepAbove) or v <= level
std::vector<ValuePoint> clipPolygon(const std::vector<ValuePoint> &poly, double level, bool keepAbove)
{
    std::vector<ValuePoint> out;
    auto inside = [&](const ValuePoint &p) { return keepAbove ? p.v >= level : p.v <= level; };
    for (std::size_t i = 0; i < poly.size(); i++)
    {
        const ValuePoint &cur = poly[i];
        const ValuePoint &prev = poly[(i + poly.size() - 1) % poly.size()];
        bool curIn = inside(cur), prevIn = inside(prev);
        if (curIn != prevIn)
        {
            double t = (level - prev.v) / (cur.v - prev.v);
            out.push_back({prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y), level});
        }
        if (curIn)
            out.push_back(cur);
    }
    return out;
}

double polygonArea(const std::vector<ValuePoint> &poly)
{
    double a = 0.0;
    for (std::size_t i = 0; i < poly.size(); i++)
    {
        const ValuePoint &p = poly[i], &q = poly[(i + 1) % poly.size()];
        a += p.x * q.y - q.x * p.y;
    }
    return 0.5 * a;
}

struct Band {
    double lo;
    double hi;
    Color color;
};

// Fills the parts of the triangles in each band (lo <= v <= hi) with the band colour
void fillBands(Figure &fig, const Matrix &coords, const std::vector<Triangle> &tris, const ColVec &values,
               const std::vector<Band> &bands)
{
    for (const auto &tri : tris)
    {
        std::vector<ValuePoint> poly;
        double vmin = inf, vmax = -inf;
        for (int n : tri)
        {
            poly.push_back({coords(n, 0), coords(n, 1), values(n)});
            vmin = std::min(vmin, values(n));
            vmax = std::max(vmax, values(n));
        }
        double triArea = std::abs(polygonArea(poly));

        for (const auto &band : bands)
        {
            if (band.hi < vmin || band.lo > vmax)
                continue;

            std::vector<ValuePoint> piece = poly;
            if (band.lo > vmin)
                piece = clipPolygon(piece, band.lo, true);
            if (piece.size() >= 3 && band.hi < vmax)
                piece = clipPolygon(piece, band.hi, false);
            if (piece.size() < 3 || std::abs(polygonArea(piece)) <= 1e-12 * triArea)
                continue;

            Polygon p;
            for (const auto &q : piece)
                p.points.push_back({q.x, q.y});
            p.fill = band.color;
            p.edge = band.color; // hides antialiasing seams between the pieces
            p.width = 1.0;
            fig.add(std::move(p));
        }
    }
}

// ---- Contour lines on triangles ----

// Contour lines of level, joined into polylines across the triangles
void addContourLines(Figure &fig, const Matrix &coords, const std::vector<Triangle> &tris, const ColVec &values,
                     double level, const Color &color, LineStyle style, double width)
{
    // Crossing points are identified by their edge (node pair)
    auto edgeKey = [&](int a, int b) {
        if (a > b)
            std::swap(a, b);
        return (long long)a * (long long)coords.rows() + b;
    };

    std::unordered_map<long long, Point> points;
    std::vector<std::array<long long, 2>> segments;

    for (const auto &tri : tris)
    {
        long long ends[2];
        int found = 0;
        for (int k = 0; k < 3; k++)
        {
            int a = tri[k], b = tri[(k + 1) % 3];
            double va = values(a), vb = values(b);
            if ((va > level) == (vb > level))
                continue;
            long long key = edgeKey(a, b);
            if (points.find(key) == points.end())
            {
                double t = (level - va) / (vb - va);
                points[key] = {coords(a, 0) + t * (coords(b, 0) - coords(a, 0)),
                               coords(a, 1) + t * (coords(b, 1) - coords(a, 1))};
            }
            if (found < 2)
                ends[found] = key;
            found++;
        }
        if (found == 2 && ends[0] != ends[1])
            segments.push_back({ends[0], ends[1]});
    }

    std::unordered_map<long long, std::vector<std::size_t>> adjacent;
    for (std::size_t s = 0; s < segments.size(); s++)
        for (long long key : segments[s])
            adjacent[key].push_back(s);

    std::vector<bool> used(segments.size(), false);
    auto next = [&](long long key) -> std::ptrdiff_t {
        for (std::size_t s : adjacent[key])
            if (!used[s])
                return std::ptrdiff_t(s);
        return -1;
    };

    for (std::size_t s = 0; s < segments.size(); s++)
    {
        if (used[s])
            continue;
        used[s] = true;
        std::deque<long long> chain = {segments[s][0], segments[s][1]};

        for (std::ptrdiff_t t; (t = next(chain.back())) >= 0;)
        {
            used[std::size_t(t)] = true;
            const auto &seg = segments[std::size_t(t)];
            chain.push_back(seg[0] == chain.back() ? seg[1] : seg[0]);
        }
        for (std::ptrdiff_t t; (t = next(chain.front())) >= 0;)
        {
            used[std::size_t(t)] = true;
            const auto &seg = segments[std::size_t(t)];
            chain.push_front(seg[0] == chain.front() ? seg[1] : seg[0]);
        }

        Line line;
        line.closed = chain.size() > 3 && chain.front() == chain.back();
        if (line.closed)
            chain.pop_back();
        for (long long key : chain)
            line.points.push_back(points[key]);
        line.color = color;
        line.style = style;
        line.width = width;
        fig.add(std::move(line));
    }
}

// ---- Shared contour drawing ----

std::vector<double> drawContourf(Figure &fig, const Matrix &coords, const std::vector<Triangle> &tris,
                                 const ColVec &values, const std::vector<double> &levels)
{
    // Band colours: the level interval midpoints (matplotlib's ContourSet layers), normalised
    // over the level range
    ColorScale scale;
    scale.colormap = fig.colormap();
    scale.vmin = levels.front();
    scale.vmax = levels.back();
    scale.boundaries = levels;

    std::vector<Band> bands;
    for (std::size_t i = 0; i + 1 < levels.size(); i++)
        bands.push_back({levels[i], levels[i + 1], scale.color(0.5 * (levels[i] + levels[i + 1]))});

    fillBands(fig, coords, tris, values, bands);
    fig.set_color_scale(scale);
    fig.set_equal_aspect(true);
    return levels;
}

std::vector<double> drawContourLines(Figure &fig, const Matrix &coords, const std::vector<Triangle> &tris,
                                     const ColVec &values, const std::vector<double> &levels,
                                     std::optional<std::pair<Color, LineStyle>> single)
{
    ColorScale scale;
    scale.colormap = fig.colormap();
    scale.vmin = levels.front();
    scale.vmax = levels.back();

    for (double level : levels)
    {
        Color color = single ? single->first : scale.color(level);
        LineStyle style = single ? single->second : LineStyle::Solid;
        addContourLines(fig, coords, tris, values, level, color, style, 1.5); // matplotlib lines.linewidth
    }
    fig.set_color_scale(scale);
    fig.set_equal_aspect(true);
    return levels;
}

// Gouraud shading approximated by 64 filled intervals between vmin and vmax; values outside get
// the end colours
void drawShaded(Figure &fig, const Matrix &coords, const std::vector<Triangle> &tris, const ColVec &values,
                double vmin, double vmax)
{
    ColorScale scale;
    scale.colormap = fig.colormap();
    scale.vmin = vmin;
    scale.vmax = vmax;

    std::vector<Band> bands;
    if (!(vmax > vmin))
        bands.push_back({-inf, inf, scale.color(vmin)});
    else
    {
        // Below / above the limits: the end colours of the map (matplotlib's under/over colours)
        bands.push_back({-inf, vmin, colormap_color(scale.colormap, 0.0)});
        const int n = 64;
        for (int i = 0; i < n; i++)
        {
            double lo = vmin + (vmax - vmin) * i / n, hi = vmin + (vmax - vmin) * (i + 1) / n;
            bands.push_back({lo, hi, scale.color(0.5 * (lo + hi))});
        }
        bands.push_back({vmax, inf, colormap_color(scale.colormap, 1.0)});
    }

    fillBands(fig, coords, tris, values, bands);
    fig.set_color_scale(scale);
    fig.set_equal_aspect(true);
}

// Mesh on top of a nodal value drawing (draw_elements in the draw_nodal_values_* functions)
void drawValueMesh(Figure &fig, const Matrix &coords, const IntMatrix &edof, int dofsPerNode, int elType,
                   const char *caller)
{
    if (dofsPerNode > 0 && elType > 0)
        addMeshLines(fig, coords, meshFaces(edof, dofsPerNode, elType, coords.rows(), caller), {0.2, 0.2, 0.2});
    else
        calfem::detail::logMessage(std::string(caller) +
                                   ": dofs_per_node and el_type must be specified to draw the mesh.");
}

std::vector<double> contourfImpl(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                 const std::vector<double> *levels, int nLevels, const std::string &title,
                                 int dofsPerNode, int elType, bool drawElements)
{
    const char *caller = "draw_nodal_values_contourf";
    checkCoords(coords, caller);
    checkNodalValues(values, coords, caller);
    auto tris = toTriangles(nodeTopology(edof, dofsPerNode, coords.rows(), caller), caller);

    auto [zmin, zmax] = valueRange(values);
    std::vector<double> lev = levels ? *levels : autoLevels(zmin, zmax, nLevels);
    checkLevels(lev, true, caller);

    drawContourf(fig, coords, tris, values, lev);
    if (drawElements)
        drawValueMesh(fig, coords, edof, dofsPerNode, elType, caller);
    setTitle(fig, title);
    return lev;
}

std::vector<double> contourImpl(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                const std::vector<double> *levels, int nLevels, const std::string &title,
                                int dofsPerNode, int elType, bool drawElements)
{
    const char *caller = "draw_nodal_values_contour";
    checkCoords(coords, caller);
    checkNodalValues(values, coords, caller);
    auto tris = toTriangles(nodeTopology(edof, dofsPerNode, coords.rows(), caller), caller);

    auto [zmin, zmax] = valueRange(values);
    std::vector<double> lev = levels ? *levels : autoLevels(zmin, zmax, nLevels);
    checkLevels(lev, false, caller);

    drawContourLines(fig, coords, tris, values, lev, std::nullopt);
    if (drawElements)
        drawValueMesh(fig, coords, edof, dofsPerNode, elType, caller);
    setTitle(fig, title);
    return lev;
}

// ---- Arrows (matplotlib's quiver) ----

// Arrow centred at (x, y) with vector (u, v), shaft width w (model units). The shape follows
// matplotlib's Quiver._h_arrows (headwidth 3, headlength 5, headaxislength 4.5, minshaft 1,
// minlength 1, pivot "middle"), in units of the shaft width.
void addArrow(Figure &fig, double x, double y, double u, double v, double w, const Color &color)
{
    const double headWidth = 3.0, headLength = 5.0, headAxisLength = 4.5, minLength = 1.0;
    const double pi = 3.14159265358979323846;

    double length = w > 0.0 ? std::hypot(u, v) / w : 0.0;
    std::vector<Point> shape;

    if (length < minLength)
    {
        // Short vectors: a hexagonal dot
        for (int k = 0; k < 6; k++)
            shape.push_back({0.5 * minLength * std::cos(k * pi / 3.0), 0.5 * minLength * std::sin(k * pi / 3.0)});
    }
    else
    {
        double minShaft = headLength; // minshaft * headlength
        double X[4], Y[4] = {0.5, 0.5, 0.5 * headWidth, 0.0};
        if (length < minShaft)
        {
            // Arrow without shaft, shrunk to the length
            double shrink = length / minShaft;
            double x0[4] = {0.0, minShaft - headAxisLength, minShaft - headLength, minShaft};
            for (int k = 0; k < 4; k++)
            {
                X[k] = shrink * x0[k];
                Y[k] *= shrink;
            }
        }
        else
        {
            X[0] = 0.0;
            X[1] = length - headAxisLength;
            X[2] = length - headLength;
            X[3] = length;
        }
        const int order[7] = {0, 1, 2, 3, 2, 1, 0};
        for (int k = 0; k < 7; k++)
            shape.push_back({X[order[k]] - 0.5 * X[3], k > 3 ? -Y[order[k]] : Y[order[k]]});
    }

    double angle = std::atan2(v, u), c = std::cos(angle), s = std::sin(angle);
    Polygon p;
    for (const auto &q : shape)
        p.points.push_back({x + w * (c * q.x - s * q.y), y + w * (s * q.x + c * q.y)});
    p.fill = color;
    p.edge = color;
    p.width = 0.5;
    fig.add(std::move(p));
}

// calfem-python's _flux_scale_factor: mean arrow length krel times the mean element size
double fluxScaleFactor(const Matrix &ex, const Matrix &ey, const Matrix &flux, double krel = 0.8)
{
    auto nel = double(ex.rows());
    double lm = 0.0, qm = 0.0;
    for (Eigen::Index e = 0; e < ex.rows(); e++)
    {
        double dx = ex.row(e).maxCoeff() - ex.row(e).minCoeff(), dy = ey.row(e).maxCoeff() - ey.row(e).minCoeff();
        lm += std::sqrt(dx * dx + dy * dy);
        qm += std::sqrt(flux(e, 0) * flux(e, 0) + flux(e, 1) * flux(e, 1));
    }
    lm /= nel;
    qm /= nel;
    return qm == 0.0 ? 0.0 : lm * krel / qm;
}

// Arrows at the element centroids; shaft width arrowWidth times the mesh width
void addFluxArrows(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &flux, double scale,
                   double arrowWidth, const Color &color, bool colorByMagnitude)
{
    double width = std::max(ex.maxCoeff() - ex.minCoeff(), ey.maxCoeff() - ey.minCoeff());
    double w = arrowWidth * width;

    ColorScale cs;
    if (colorByMagnitude)
    {
        Eigen::VectorXd magnitude = flux.leftCols(2).rowwise().norm();
        cs.colormap = fig.colormap();
        cs.vmin = magnitude.minCoeff();
        cs.vmax = magnitude.maxCoeff();
    }

    for (Eigen::Index e = 0; e < ex.rows(); e++)
    {
        Color c = colorByMagnitude ? cs.color(std::hypot(flux(e, 0), flux(e, 1))) : color;
        addArrow(fig, ex.row(e).mean(), ey.row(e).mean(), scale * flux(e, 0), scale * flux(e, 1), w, c);
    }

    if (colorByMagnitude)
        fig.set_color_scale(cs);
    fig.set_equal_aspect(true);
}

// Element coordinates of the faces
std::pair<Matrix, Matrix> faceCoords(const Matrix &coords, const std::vector<std::vector<int>> &faces)
{
    auto nen = faces.empty() ? 0 : Eigen::Index(faces[0].size());
    Matrix ex(Eigen::Index(faces.size()), nen), ey(Eigen::Index(faces.size()), nen);
    for (std::size_t e = 0; e < faces.size(); e++)
        for (Eigen::Index k = 0; k < nen; k++)
        {
            ex(Eigen::Index(e), k) = coords(faces[e][std::size_t(k)], 0);
            ey(Eigen::Index(e), k) = coords(faces[e][std::size_t(k)], 1);
        }
    return {ex, ey};
}

// ---- eliso2 ----

std::vector<double> eliso2Impl(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed,
                               const std::vector<double> *isov, int nIsov,
                               std::optional<std::array<int, 2>> plotpar)
{
    const char *caller = "eliso2";
    if (ex.rows() != ey.rows() || ex.cols() != ey.cols())
        throw std::invalid_argument("eliso2: check size of ex, ey dimensions.");
    if (ed.rows() != ex.rows() || ed.cols() != ex.cols())
        throw std::invalid_argument("eliso2: ed must have one value per element node (same size as ex).");

    std::optional<std::pair<Color, LineStyle>> single;
    if (plotpar)
        single = std::make_pair(codes::colorCode((*plotpar)[1], "plotpar[1]"),
                                codes::lineCode((*plotpar)[0], "plotpar[0]"));

    // Merge the element nodes into a node based mesh (as calfem-python, with a relative tolerance)
    double xr = ex.size() ? ex.maxCoeff() - ex.minCoeff() : 0.0;
    double yr = ey.size() ? ey.maxCoeff() - ey.minCoeff() : 0.0;
    double tol = 1e-9 * std::max(xr, yr);
    if (!(tol > 0.0))
        tol = 1.0;

    std::map<std::pair<long long, long long>, int> index;
    std::vector<Point> nodes;
    std::vector<double> nodeValues;
    std::vector<std::vector<int>> topo(static_cast<std::size_t>(ex.rows()));
    for (Eigen::Index e = 0; e < ex.rows(); e++)
        for (Eigen::Index k = 0; k < ex.cols(); k++)
        {
            std::pair<long long, long long> key{std::llround(ex(e, k) / tol), std::llround(ey(e, k) / tol)};
            auto it = index.find(key);
            int n;
            if (it == index.end())
            {
                n = int(nodes.size());
                index[key] = n;
                nodes.push_back({ex(e, k), ey(e, k)});
                nodeValues.push_back(ed(e, k));
            }
            else
                n = it->second;
            topo[static_cast<std::size_t>(e)].push_back(n);
        }

    Matrix coords(Eigen::Index(nodes.size()), 2);
    ColVec values(Eigen::Index(nodes.size()));
    for (std::size_t i = 0; i < nodes.size(); i++)
    {
        coords(Eigen::Index(i), 0) = nodes[i].x;
        coords(Eigen::Index(i), 1) = nodes[i].y;
        values(Eigen::Index(i)) = nodeValues[i];
    }

    auto tris = toTriangles(topo, caller);
    auto [zmin, zmax] = valueRange(values);
    std::vector<double> lev = isov ? *isov : autoLevels(zmin, zmax, nIsov);
    checkLevels(lev, false, caller);
    return drawContourLines(fig, coords, tris, values, lev, single);
}

} // namespace

// ---- draw_mesh ----

void draw_mesh(const Matrix &coords, const IntMatrix &edof, int dofs_per_node, int el_type, const std::string &title,
               const Color &color, const Color &face_color, bool filled, bool show_nodes, const Color &node_color)
{
    draw_mesh(gcf(), coords, edof, dofs_per_node, el_type, title, color, face_color, filled, show_nodes, node_color);
}

void draw_mesh(Figure &fig, const Matrix &coords, const IntMatrix &edof, int dofs_per_node, int el_type,
               const std::string &title, const Color &color, const Color &face_color, bool filled, bool show_nodes,
               const Color &node_color)
{
    checkCoords(coords, "draw_mesh");
    auto faces = meshFaces(edof, dofs_per_node, el_type, coords.rows(), "draw_mesh");

    if (filled)
    {
        for (const auto &face : faces)
        {
            Polygon p;
            p.points = facePoints(coords, face);
            p.fill = face_color;
            p.edge = color;
            fig.add(std::move(p));
        }
    }
    else
        addMeshLines(fig, coords, faces, color);

    if (show_nodes)
    {
        Markers m;
        for (Eigen::Index i = 0; i < coords.rows(); i++)
            m.points.push_back({coords(i, 0), coords(i, 1)});
        m.marker = Marker::Circle;
        m.color = node_color;
        m.filled = true;
        fig.add(std::move(m));
    }

    fig.set_equal_aspect(true);
    setTitle(fig, title);
}

// ---- draw_elements ----

void draw_elements(const Matrix &ex, const Matrix &ey, const std::string &title, const Color &color,
                   const Color &face_color, const Color &node_color, LineStyle line_style, bool filled, bool closed,
                   bool show_nodes)
{
    draw_elements(gcf(), ex, ey, title, color, face_color, node_color, line_style, filled, closed, show_nodes);
}

void draw_elements(Figure &fig, const Matrix &ex, const Matrix &ey, const std::string &title, const Color &color,
                   const Color &face_color, const Color &node_color, LineStyle line_style, bool filled, bool closed,
                   bool show_nodes)
{
    if (ex.rows() != ey.rows() || ex.cols() != ey.cols())
        throw std::invalid_argument("draw_elements: check size of ex, ey dimensions.");

    for (Eigen::Index e = 0; e < ex.rows(); e++)
    {
        std::vector<Point> points;
        for (Eigen::Index k = 0; k < ex.cols(); k++)
            points.push_back({ex(e, k), ey(e, k)});

        if (filled)
        {
            Polygon p;
            p.points = std::move(points);
            p.fill = face_color;
            p.edge = color;
            p.style = line_style;
            fig.add(std::move(p));
        }
        else
        {
            Line l;
            l.points = std::move(points);
            l.color = color;
            l.style = line_style;
            l.closed = closed;
            fig.add(std::move(l));
        }
    }

    if (show_nodes)
        draw_node_circles(fig, ex, ey, "", node_color, face_color, true);

    fig.set_equal_aspect(true);
    setTitle(fig, title);
}

// ---- draw_node_circles ----

void draw_node_circles(const Matrix &ex, const Matrix &ey, const std::string &title, const Color &color,
                       const Color &face_color, bool filled, Marker marker_type)
{
    draw_node_circles(gcf(), ex, ey, title, color, face_color, filled, marker_type);
}

void draw_node_circles(Figure &fig, const Matrix &ex, const Matrix &ey, const std::string &title, const Color &color,
                       const Color & /*face_color*/, bool filled, Marker marker_type)
{
    if (ex.rows() != ey.rows() || ex.cols() != ey.cols())
        throw std::invalid_argument("draw_node_circles: check size of ex, ey dimensions.");

    // As calfem-python: filled markers in color, otherwise outlined (face_color is not used)
    Markers m;
    for (Eigen::Index e = 0; e < ex.rows(); e++)
        for (Eigen::Index k = 0; k < ex.cols(); k++)
            m.points.push_back({ex(e, k), ey(e, k)});
    m.marker = marker_type;
    m.color = color;
    m.filled = filled;
    fig.add(std::move(m));

    fig.set_equal_aspect(true);
    setTitle(fig, title);
}

// ---- draw_element_values ----

void draw_element_values(const ColVec &values, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                         int el_type, const Matrix &displacements, bool draw_elements, bool draw_undisplaced_mesh,
                         double magnfac, const std::string &title, const Color &color, const Color &node_color)
{
    draw_element_values(gcf(), values, coords, edof, dofs_per_node, el_type, displacements, draw_elements,
                        draw_undisplaced_mesh, magnfac, title, color, node_color);
}

void draw_element_values(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                         int dofs_per_node, int el_type, const Matrix &displacements, bool draw_elements,
                         bool draw_undisplaced_mesh, double magnfac, const std::string &title, const Color &color,
                         const Color & /*node_color*/)
{
    const char *caller = "draw_element_values";
    checkCoords(coords, caller);
    auto faces = meshFaces(edof, dofs_per_node, el_type, coords.rows(), caller);
    if (values.size() != Eigen::Index(faces.size()))
        throw std::invalid_argument("draw_element_values: there must be one value per element.");

    if (draw_undisplaced_mesh)
        addMeshLines(fig, coords, faces, {0.5, 0.5, 0.5});

    Matrix xy = coords;
    if (displacements.size() > 0)
        xy += magnfac * nodalDisplacements(displacements, coords, dofs_per_node, caller);

    ColorScale scale;
    scale.colormap = fig.colormap();
    auto [vmin, vmax] = valueRange(values);
    scale.vmin = vmin;
    scale.vmax = vmax;

    for (std::size_t e = 0; e < faces.size(); e++)
    {
        Polygon p;
        p.points = facePoints(xy, faces[e]);
        p.fill = scale.color(values(Eigen::Index(e)));
        if (draw_elements)
            p.edge = color;
        fig.add(std::move(p));
    }

    fig.set_color_scale(scale);
    fig.set_equal_aspect(true);
    setTitle(fig, title);
}

// ---- draw_displacements ----

double draw_displacements(const Matrix &a, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                          int el_type, bool draw_undisplaced_mesh, double magnfac, double magscale,
                          const std::string &title, const Color &color)
{
    return draw_displacements(gcf(), a, coords, edof, dofs_per_node, el_type, draw_undisplaced_mesh, magnfac,
                              magscale, title, color);
}

double draw_displacements(Figure &fig, const Matrix &a, const Matrix &coords, const IntMatrix &edof,
                          int dofs_per_node, int el_type, bool draw_undisplaced_mesh, double magnfac, double magscale,
                          const std::string &title, const Color &color)
{
    const char *caller = "draw_displacements";
    checkCoords(coords, caller);
    auto faces = meshFaces(edof, dofs_per_node, el_type, coords.rows(), caller);

    if (draw_undisplaced_mesh)
        addMeshLines(fig, coords, faces, {0.8, 0.8, 0.8});

    Matrix xy = coords;
    if (a.size() > 0)
    {
        Matrix u = nodalDisplacements(a, coords, dofs_per_node, caller);
        if (magnfac < 0.0)
        {
            double maxSize = (coords.colwise().maxCoeff() - coords.colwise().minCoeff()).maxCoeff();
            double maxDisp = u.rowwise().norm().maxCoeff();
            magnfac = maxDisp > 0.0 ? magscale * maxSize / maxDisp : 1.0;
        }
        xy += magnfac * u;
    }
    else if (magnfac < 0.0)
        magnfac = 1.0;

    addMeshLines(fig, xy, faces, color);
    fig.set_equal_aspect(true);
    setTitle(fig, title);
    return magnfac;
}

// ---- draw_element_flux, elflux2 ----

double draw_element_flux(const Matrix &flux, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                         int el_type, std::optional<double> scale, const Color &color, bool color_by_magnitude,
                         bool draw_elements, const Color &element_color, double arrow_width, const std::string &title)
{
    return draw_element_flux(gcf(), flux, coords, edof, dofs_per_node, el_type, scale, color, color_by_magnitude,
                             draw_elements, element_color, arrow_width, title);
}

double draw_element_flux(Figure &fig, const Matrix &flux, const Matrix &coords, const IntMatrix &edof,
                         int dofs_per_node, int el_type, std::optional<double> scale, const Color &color,
                         bool color_by_magnitude, bool draw_elements, const Color &element_color,
                         double arrow_width, const std::string &title)
{
    const char *caller = "draw_element_flux";
    checkCoords(coords, caller);
    if (flux.cols() < 2)
        throw std::invalid_argument("draw_element_flux: flux must be an E x 2 array of [qx, qy].");
    if (flux.rows() != edof.rows())
        throw std::invalid_argument("draw_element_flux: there must be one row of flux for each element.");

    auto faces = meshFaces(edof, dofs_per_node, el_type, coords.rows(), caller);
    auto [ex, ey] = faceCoords(coords, faces);

    double s = scale ? *scale : fluxScaleFactor(ex, ey, flux);

    if (draw_elements)
        addMeshLines(fig, coords, faces, element_color);
    addFluxArrows(fig, ex, ey, flux, s, arrow_width, color, color_by_magnitude);
    setTitle(fig, title);
    return s;
}

double elflux2(const Matrix &ex, const Matrix &ey, const Matrix &es, int plotcolor, std::optional<double> sfac)
{
    return elflux2(gcf(), ex, ey, es, plotcolor, sfac);
}

double elflux2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &es, int plotcolor,
               std::optional<double> sfac)
{
    if (ex.rows() != ey.rows() || ex.cols() != ey.cols())
        throw std::invalid_argument("elflux2: check size of coordinate input arguments.");
    if (es.rows() != ex.rows())
        throw std::invalid_argument("elflux2: there must be one row of es for each element.");
    if (ex.cols() != 3 && ex.cols() != 4)
        throw std::invalid_argument("elflux2: only 3- and 4-node elements are supported.");
    if (es.cols() < 2)
        throw std::invalid_argument("elflux2: es must contain the x- and y-components of the flow.");

    double s = sfac ? *sfac : fluxScaleFactor(ex, ey, es);

    // As calfem-python, an unknown colour code gives black
    Color color = plotcolor >= 1 && plotcolor <= 4 ? codes::colorCode(plotcolor, "plotcolor") : Color{0.0, 0.0, 0.0};
    addFluxArrows(fig, ex, ey, es, s, 0.003, color, false);
    return s;
}

// ---- draw_nodal_values_contourf, _contour, _shaded, draw_nodal_values ----

std::vector<double> draw_nodal_values_contourf(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                               int levels, const std::string &title, int dofs_per_node, int el_type,
                                               bool draw_elements)
{
    return contourfImpl(gcf(), values, coords, edof, nullptr, levels, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contourf(Figure &fig, const ColVec &values, const Matrix &coords,
                                               const IntMatrix &edof, int levels, const std::string &title,
                                               int dofs_per_node, int el_type, bool draw_elements)
{
    return contourfImpl(fig, values, coords, edof, nullptr, levels, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contourf(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                               const std::vector<double> &levels, const std::string &title,
                                               int dofs_per_node, int el_type, bool draw_elements)
{
    return contourfImpl(gcf(), values, coords, edof, &levels, 0, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contourf(Figure &fig, const ColVec &values, const Matrix &coords,
                                               const IntMatrix &edof, const std::vector<double> &levels,
                                               const std::string &title, int dofs_per_node, int el_type,
                                               bool draw_elements)
{
    return contourfImpl(fig, values, coords, edof, &levels, 0, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contour(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                              int levels, const std::string &title, int dofs_per_node, int el_type,
                                              bool draw_elements)
{
    return contourImpl(gcf(), values, coords, edof, nullptr, levels, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contour(Figure &fig, const ColVec &values, const Matrix &coords,
                                              const IntMatrix &edof, int levels, const std::string &title,
                                              int dofs_per_node, int el_type, bool draw_elements)
{
    return contourImpl(fig, values, coords, edof, nullptr, levels, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contour(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                              const std::vector<double> &levels, const std::string &title,
                                              int dofs_per_node, int el_type, bool draw_elements)
{
    return contourImpl(gcf(), values, coords, edof, &levels, 0, title, dofs_per_node, el_type, draw_elements);
}

std::vector<double> draw_nodal_values_contour(Figure &fig, const ColVec &values, const Matrix &coords,
                                              const IntMatrix &edof, const std::vector<double> &levels,
                                              const std::string &title, int dofs_per_node, int el_type,
                                              bool draw_elements)
{
    return contourImpl(fig, values, coords, edof, &levels, 0, title, dofs_per_node, el_type, draw_elements);
}

void draw_nodal_values_shaded(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                              const std::string &title, int dofs_per_node, int el_type, bool draw_elements)
{
    draw_nodal_values_shaded(gcf(), values, coords, edof, title, dofs_per_node, el_type, draw_elements);
}

void draw_nodal_values_shaded(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                              const std::string &title, int dofs_per_node, int el_type, bool draw_elements)
{
    const char *caller = "draw_nodal_values_shaded";
    checkCoords(coords, caller);
    checkNodalValues(values, coords, caller);
    auto tris = toTriangles(nodeTopology(edof, dofs_per_node, coords.rows(), caller), caller);

    auto [vmin, vmax] = valueRange(values);
    drawShaded(fig, coords, tris, values, vmin, vmax);
    if (draw_elements)
        drawValueMesh(fig, coords, edof, dofs_per_node, el_type, caller);
    setTitle(fig, title);
}

void draw_nodal_values(const ColVec &values, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                       int el_type, std::optional<std::array<double, 2>> clim, bool draw_elements,
                       const std::string &title)
{
    draw_nodal_values(gcf(), values, coords, edof, dofs_per_node, el_type, clim, draw_elements, title);
}

void draw_nodal_values(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                       int dofs_per_node, int el_type, std::optional<std::array<double, 2>> clim, bool draw_elements,
                       const std::string &title)
{
    const char *caller = "draw_nodal_values";
    checkCoords(coords, caller);
    checkNodalValues(values, coords, caller);
    if (dofs_per_node < 1)
        throw std::invalid_argument("draw_nodal_values: dofs_per_node must be positive.");
    auto tris = toTriangles(nodeTopology(edof, dofs_per_node, coords.rows(), caller), caller);

    auto [vmin, vmax] = valueRange(values);
    if (clim)
    {
        vmin = (*clim)[0];
        vmax = (*clim)[1];
    }
    drawShaded(fig, coords, tris, values, vmin, vmax);
    if (draw_elements)
        addMeshLines(fig, coords, meshFaces(edof, dofs_per_node, el_type, coords.rows(), caller), {0.2, 0.2, 0.2});
    setTitle(fig, title);
}

// ---- eliso2 ----

std::vector<double> eliso2(const Matrix &ex, const Matrix &ey, const Matrix &ed, int isov,
                           std::optional<std::array<int, 2>> plotpar)
{
    return eliso2Impl(gcf(), ex, ey, ed, nullptr, isov, plotpar);
}

std::vector<double> eliso2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed, int isov,
                           std::optional<std::array<int, 2>> plotpar)
{
    return eliso2Impl(fig, ex, ey, ed, nullptr, isov, plotpar);
}

std::vector<double> eliso2(const Matrix &ex, const Matrix &ey, const Matrix &ed, const std::vector<double> &isov,
                           std::optional<std::array<int, 2>> plotpar)
{
    return eliso2Impl(gcf(), ex, ey, ed, &isov, 0, plotpar);
}

std::vector<double> eliso2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed,
                           const std::vector<double> &isov, std::optional<std::array<int, 2>> plotpar)
{
    return eliso2Impl(fig, ex, ey, ed, &isov, 0, plotpar);
}

// ---- elval2 ----

void elval2(const Matrix &ex, const Matrix &ey, const ColVec &ev, bool showMesh)
{
    elval2(gcf(), ex, ey, ev, showMesh);
}

void elval2(Figure &fig, const Matrix &ex, const Matrix &ey, const ColVec &ev, bool showMesh)
{
    if (ex.rows() != ey.rows() || ex.cols() != ey.cols())
        throw std::invalid_argument("elval2: check size of ex, ey dimensions.");
    if (ev.size() != ex.rows())
        throw std::invalid_argument("elval2: there must be one value per element.");

    ColorScale scale;
    scale.colormap = fig.colormap();
    auto [vmin, vmax] = valueRange(ev);
    scale.vmin = vmin;
    scale.vmax = vmax;

    for (Eigen::Index e = 0; e < ex.rows(); e++)
    {
        Polygon p;
        for (Eigen::Index k = 0; k < ex.cols(); k++)
            p.points.push_back({ex(e, k), ey(e, k)});
        p.fill = scale.color(ev(e));
        p.edge = showMesh ? Color{0.0, 0.0, 0.0} : p.fill; // matplotlib edgecolor "face" without the mesh
        p.width = 0.5;
        fig.add(std::move(p));
    }

    fig.set_color_scale(scale);
    fig.set_equal_aspect(true);
}

} // namespace calfem::vis
