//
// Tests of the calfem::vis colour maps, colour bar and mesh/field drawing functions
// (comparisons with calfem-python are in test_vis_fields_reference.cpp).
//

#include <calfem/vis.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace calfem;
using Catch::Approx;

namespace {

int count(const std::string &text, const std::string &pattern)
{
    int n = 0;
    for (auto pos = text.find(pattern); pos != std::string::npos; pos = text.find(pattern, pos + 1))
        n++;
    return n;
}

// nx x ny quadrilaterals on [0, 2] x [0, 1]
struct Grid {
    Matrix coords;
    IntMatrix topo; // 1-based node numbers

    // Element dofs with dofsPerNode dofs per node
    IntMatrix edof(int dofsPerNode) const
    {
        IntMatrix e(topo.rows(), topo.cols() * dofsPerNode);
        for (Eigen::Index i = 0; i < topo.rows(); i++)
            for (Eigen::Index k = 0; k < topo.cols(); k++)
                for (int d = 0; d < dofsPerNode; d++)
                    e(i, k * dofsPerNode + d) = (topo(i, k) - 1) * dofsPerNode + d + 1;
        return e;
    }

    // Field f(x, y) at the nodes
    template <typename F> ColVec field(F f) const
    {
        ColVec v(coords.rows());
        for (Eigen::Index i = 0; i < coords.rows(); i++)
            v(i) = f(coords(i, 0), coords(i, 1));
        return v;
    }
};

Grid grid(int nx, int ny)
{
    Grid g;
    g.coords.resize((nx + 1) * (ny + 1), 2);
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
            g.coords.row(j * (nx + 1) + i) << 2.0 * i / nx, 1.0 * j / ny;
    g.topo.resize(nx * ny, 4);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
            int n = j * (nx + 1) + i + 1;
            g.topo.row(j * nx + i) << n, n + 1, n + nx + 2, n + nx + 1;
        }
    return g;
}

double area(const std::vector<vis::Point> &points)
{
    double a = 0.0;
    for (std::size_t i = 0; i < points.size(); i++)
    {
        const auto &p = points[i], &q = points[(i + 1) % points.size()];
        a += p.x * q.y - q.x * p.y;
    }
    return 0.5 * std::abs(a);
}

std::vector<vis::Polygon> polygons(const vis::Figure &fig)
{
    std::vector<vis::Polygon> result;
    for (const auto &item : fig.items())
        if (auto p = std::get_if<vis::Polygon>(&item))
            result.push_back(*p);
    return result;
}

std::vector<vis::Line> lines(const vis::Figure &fig)
{
    std::vector<vis::Line> result;
    for (const auto &item : fig.items())
        if (auto l = std::get_if<vis::Line>(&item))
            result.push_back(*l);
    return result;
}

bool sameColor(const vis::Color &a, const vis::Color &b)
{
    return std::abs(a.r - b.r) < 1e-12 && std::abs(a.g - b.g) < 1e-12 && std::abs(a.b - b.b) < 1e-12;
}

double totalArea(const vis::Figure &fig, const vis::Color &color)
{
    double a = 0.0;
    for (const auto &p : polygons(fig))
        if (sameColor(p.fill, color))
            a += area(p.points);
    return a;
}

} // namespace

TEST_CASE("Colour maps are matplotlib's lookup tables", "[vis][fields]")
{
    auto check = [](const vis::Color &c, double r, double g, double b) {
        CHECK(c.r == Approx(r));
        CHECK(c.g == Approx(g));
        CHECK(c.b == Approx(b));
    };

    check(vis::colormap_color(vis::Colormap::Viridis, 0.0), 0.267004, 0.004874, 0.329415);
    check(vis::colormap_color(vis::Colormap::Viridis, 0.5), 0.127568, 0.566949, 0.550556); // entry 128
    check(vis::colormap_color(vis::Colormap::Viridis, 1.0), 0.993248, 0.906157, 0.143936);
    check(vis::colormap_color(vis::Colormap::Viridis, -0.5), 0.267004, 0.004874, 0.329415); // clamped
    check(vis::colormap_color(vis::Colormap::Viridis, 2.0), 0.993248, 0.906157, 0.143936);

    check(vis::colormap_color(vis::Colormap::Jet, 0.0), 0.0, 0.0, 0.5);
    check(vis::colormap_color(vis::Colormap::Jet, 1.0), 0.5, 0.0, 0.0);
    check(vis::colormap_color(vis::Colormap::Jet, 0.5), 0.4901960784313725, 1.0, 0.4775458570524984);
    check(vis::colormap_color(vis::Colormap::Gray, 1.0), 1.0, 1.0, 1.0);

    // An empty value range maps everything to the first colour (matplotlib's Normalize)
    vis::ColorScale scale;
    scale.vmin = scale.vmax = 3.0;
    check(scale.color(5.0), 0.267004, 0.004874, 0.329415);
}

TEST_CASE("Filled contours cover the level intervals", "[vis][fields]")
{
    Grid g = grid(4, 2);
    ColVec v = g.field([](double x, double) { return x; });

    vis::Figure fig;
    auto levels = vis::draw_nodal_values_contourf(fig, v, g.coords, g.topo, std::vector<double>{0.0, 0.5, 1.0, 1.5, 2.0});
    REQUIRE(levels.size() == 5);

    // f = x: every interval is a 0.5 x 1 strip, coloured by its midpoint over the level range
    const auto &scale = fig.color_scale();
    REQUIRE(scale);
    CHECK(scale->boundaries == levels);
    for (int i = 0; i < 4; i++)
    {
        vis::Color c = vis::colormap_color(vis::Colormap::Viridis, (0.25 + 0.5 * i) / 2.0);
        CHECK(totalArea(fig, c) == Approx(0.5));
    }

    // Regions outside the levels are not filled
    vis::Figure partial;
    vis::draw_nodal_values_contourf(partial, v, g.coords, g.topo, std::vector<double>{0.5, 1.0});
    double a = 0.0;
    for (const auto &p : polygons(partial))
        a += area(p.points);
    CHECK(a == Approx(0.5));
}

TEST_CASE("Automatic contour levels follow matplotlib", "[vis][fields]")
{
    Grid g = grid(2, 2);
    vis::Figure fig;
    CHECK(vis::draw_nodal_values_contourf(fig, g.field([](double x, double y) { return x + y; }), g.coords, g.topo,
                                          5) == std::vector<double>{0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0});

    // A constant field still gives a level interval around the value
    auto levels = vis::draw_nodal_values_contourf(fig, g.field([](double, double) { return 2.0; }), g.coords, g.topo);
    REQUIRE(levels.size() >= 2);
    CHECK(levels.front() <= 2.0);
    CHECK(levels.back() >= 2.0);
}

TEST_CASE("Contour lines are joined across the elements", "[vis][fields]")
{
    Grid g = grid(4, 4);

    // f = x: the level 0.75 is a vertical line through all elements
    vis::Figure fig;
    vis::draw_nodal_values_contour(fig, g.field([](double x, double) { return x; }), g.coords, g.topo,
                                   std::vector<double>{0.75});
    auto l = lines(fig);
    REQUIRE(l.size() == 1);
    CHECK_FALSE(l[0].closed);
    CHECK(l[0].points.size() == 2 * 4 + 1); // two triangles per element row, plus the end point
    for (const auto &p : l[0].points)
        CHECK(p.x == Approx(0.75));

    // A circle around the centre gives one closed line
    vis::Figure circle;
    vis::draw_nodal_values_contour(circle, g.field([](double x, double y) { return std::hypot(x - 1.0, y - 0.5); }),
                                   g.coords, g.topo, std::vector<double>{0.3});
    auto c = lines(circle);
    REQUIRE(c.size() == 1);
    CHECK(c[0].closed);
}

TEST_CASE("Element dofs are converted to node topology", "[vis][fields]")
{
    Grid g = grid(3, 2);
    ColVec v = g.field([](double x, double y) { return x * x + y; });

    vis::Figure byTopology, byDofs;
    vis::draw_nodal_values_contourf(byTopology, v, g.coords, g.topo, 6);
    vis::draw_nodal_values_contourf(byDofs, v, g.coords, g.edof(2), 6, "", 2, 3);

    auto a = polygons(byTopology), b = polygons(byDofs);
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); i++)
        CHECK(area(a[i].points) == Approx(area(b[i].points)));
}

TEST_CASE("Shaded values use the colour limits", "[vis][fields]")
{
    Grid g = grid(4, 2);
    ColVec v = g.field([](double x, double) { return x; });

    vis::Figure fig;
    vis::draw_nodal_values(fig, v, g.coords, g.edof(1), 1, 3, std::array<double, 2>{0.5, 1.5}, false);
    REQUIRE(fig.color_scale());
    CHECK(fig.color_scale()->vmin == 0.5);
    CHECK(fig.color_scale()->vmax == 1.5);

    // Values below / above the limits get the end colours; the whole mesh is filled
    CHECK(totalArea(fig, vis::colormap_color(vis::Colormap::Viridis, 0.0)) >= 0.5 - 1e-12);
    CHECK(totalArea(fig, vis::colormap_color(vis::Colormap::Viridis, 1.0)) >= 0.5 - 1e-12);
    double a = 0.0;
    for (const auto &p : polygons(fig))
        a += area(p.points);
    CHECK(a == Approx(2.0));
}

TEST_CASE("The figure colour map is used", "[vis][fields]")
{
    Grid g = grid(2, 1);
    vis::Figure fig;
    vis::colormap(fig, vis::Colormap::Jet);
    vis::draw_element_values(fig, ColVec::LinSpaced(2, 0.0, 1.0), g.coords, g.edof(1), 1, 3);

    auto p = polygons(fig);
    REQUIRE(p.size() == 2);
    CHECK(sameColor(p[0].fill, vis::colormap_color(vis::Colormap::Jet, 0.0)));
    CHECK(sameColor(p[1].fill, vis::colormap_color(vis::Colormap::Jet, 1.0)));
    REQUIRE(p[0].edge);
    CHECK(sameColor(*p[0].edge, {0.0, 0.0, 0.0}));
}

TEST_CASE("Colour bar layout and output", "[vis][fields]")
{
    Grid g = grid(4, 2);
    ColVec v = g.field([](double x, double y) { return 10.0 * x + y; });

    vis::Figure empty;
    CHECK_THROWS_AS(vis::colorbar(empty), std::runtime_error);

    vis::Figure fig;
    vis::draw_nodal_values_shaded(fig, v, g.coords, g.topo);
    vis::Layout without = vis::layout(fig, 800.0, 600.0);
    CHECK_FALSE(without.colorbar);

    vis::colorbar(fig);
    REQUIRE(fig.colorbar().has_value());
    CHECK(fig.colorbar()->vmax == fig.color_scale()->vmax);
    vis::Layout l = vis::layout(fig, 800.0, 600.0);
    REQUIRE(l.colorbar);
    CHECK(l.plotWidth < without.plotWidth);
    CHECK(l.colorbarLeft > l.plotLeft + l.plotWidth);
    CHECK(l.colorbarLeft + l.colorbarWidth < 800.0);
    CHECK(l.colorbarMin == v.minCoeff());
    CHECK(l.colorbarMax == v.maxCoeff());
    CHECK(l.cy(l.colorbarMax) == Approx(l.colorbarTop));
    CHECK(l.cy(l.colorbarMin) == Approx(l.colorbarTop + l.colorbarHeight));
    REQUIRE_FALSE(l.colorbarTicks.empty());
    for (double t : l.colorbarTicks)
        CHECK((t >= l.colorbarMin && t <= l.colorbarMax));

    std::string svg = vis::to_svg(fig, 800.0, 600.0);
    CHECK(count(svg, "<polygon") == int(polygons(fig).size()));
    CHECK(count(svg, "<rect") >= 128); // background, colour bar strips and frame

    // Filled contours: one band per level interval, ticks at the levels
    vis::Figure contours;
    auto levels = vis::draw_nodal_values_contourf(contours, v, g.coords, g.topo, 4);
    vis::colorbar(contours);
    vis::Layout lc = vis::layout(contours, 800.0, 600.0);
    CHECK(lc.colorbarTicks == levels);
    CHECK(lc.colorbarMin == levels.front());
    CHECK(lc.colorbarMax == levels.back());

    // The colour bar keeps the scale of the drawing it was made for
    vis::draw_nodal_values_contour(contours, v, g.coords, g.topo, 3);
    CHECK(contours.colorbar()->boundaries == levels);
    CHECK(contours.color_scale()->boundaries.empty());

    // clear() removes the colour scale and the colour bar
    contours.clear();
    CHECK_FALSE(contours.color_scale());
    CHECK_FALSE(contours.colorbar().has_value());
}

TEST_CASE("Mesh, displacement and flux drawings", "[vis][fields]")
{
    Grid g = grid(2, 1);

    vis::Figure mesh;
    vis::draw_mesh(mesh, g.coords, g.edof(2), 2, 3, "Mesh", {0, 0, 1}, {0.8, 0.8, 0.8}, false, true);
    CHECK(lines(mesh).size() == 2);
    CHECK(mesh.title() == "Mesh");

    // Automatic magnification: the largest displacement becomes 0.1 times the model size
    ColVec a = ColVec::Zero(2 * g.coords.rows());
    a(2 * 5 + 1) = 0.01; // uy of the last node
    vis::Figure disp;
    double magnfac = vis::draw_displacements(disp, a, g.coords, g.edof(2), 2, 3);
    CHECK(magnfac == Approx(0.1 * 2.0 / 0.01));
    auto l = lines(disp);
    REQUIRE(l.size() == 2);
    CHECK(l[1].points[2].y == Approx(1.0 + 0.2));

    // Flux arrows: one polygon per element, default scale from the mean element size
    Matrix flux(2, 2);
    flux << 1.0, 0.0, 0.0, 2.0;
    vis::Figure arrows;
    double scale = vis::draw_element_flux(arrows, flux, g.coords, g.edof(1), 1, 3);
    CHECK(scale == Approx(0.8 * std::sqrt(2.0) / 1.5));
    CHECK(polygons(arrows).size() == 2);

    // Zero flux gives a zero scale (calfem-python)
    vis::Figure none;
    CHECK(vis::draw_element_flux(none, Matrix::Zero(2, 2), g.coords, g.edof(1), 1, 3) == 0.0);
}

TEST_CASE("eliso2 merges the element nodes", "[vis][fields]")
{
    // Two quadrilaterals sharing an edge; f = x gives one line through both
    Matrix ex(2, 4), ey(2, 4), ed(2, 4);
    ex << 0, 1, 1, 0, 1, 2, 2, 1;
    ey << 0, 0, 1, 1, 0, 0, 1, 1;
    ed = ex;

    vis::Figure fig;
    vis::eliso2(fig, ex, ey, ed, std::vector<double>{0.5, 1.5}, std::array<int, 2>{2, 4});
    auto l = lines(fig);
    REQUIRE(l.size() == 2);
    for (const auto &line : l)
    {
        CHECK(line.style == vis::LineStyle::Dashed);
        CHECK(sameColor(line.color, {1.0, 0.0, 0.0}));
    }

    // f = y: the level 0.5 runs through both elements as one line
    vis::Figure across;
    vis::eliso2(across, ex, ey, ey, std::vector<double>{0.5});
    auto a = lines(across);
    REQUIRE(a.size() == 1);
    double x0 = a[0].points.front().x, x1 = a[0].points.back().x;
    CHECK(std::min(x0, x1) == Approx(0.0).margin(1e-12));
    CHECK(std::max(x0, x1) == Approx(2.0));
}

TEST_CASE("Invalid field drawing inputs throw", "[vis][fields]")
{
    Grid g = grid(2, 1);
    ColVec v = g.field([](double x, double) { return x; });
    vis::Figure fig;

    CHECK_THROWS_AS(vis::draw_nodal_values_contourf(fig, ColVec::Zero(2), g.coords, g.topo), std::invalid_argument);
    CHECK_THROWS_AS(vis::draw_nodal_values_contourf(fig, v, g.coords, g.topo, std::vector<double>{1.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(vis::draw_nodal_values_contour(fig, v, g.coords, g.topo, std::vector<double>{1.0, 0.5}),
                    std::invalid_argument);
    CHECK_THROWS_AS(vis::draw_mesh(fig, g.coords, g.edof(1), 1, 5), std::invalid_argument);
    CHECK_THROWS_AS(vis::draw_element_values(fig, ColVec::Zero(3), g.coords, g.edof(1), 1, 3), std::invalid_argument);
    CHECK_THROWS_AS(vis::eliso2(fig, Matrix::Zero(1, 4), Matrix::Zero(1, 4), Matrix::Zero(1, 3)),
                    std::invalid_argument);
    CHECK_THROWS_AS(vis::elflux2(fig, Matrix::Zero(1, 5), Matrix::Zero(1, 5), Matrix::Zero(1, 2)),
                    std::invalid_argument);
    CHECK_THROWS_AS(vis::elval2(fig, Matrix::Zero(2, 4), Matrix::Zero(2, 4), ColVec::Zero(3)), std::invalid_argument);

    IntMatrix outside = g.topo;
    outside(0, 0) = 99;
    CHECK_THROWS_AS(vis::draw_nodal_values_shaded(fig, v, g.coords, outside), std::invalid_argument);
}
