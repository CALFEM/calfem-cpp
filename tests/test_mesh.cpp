//
// Geometry and mesh generation (built-in mesher).
//

#include <calfem/calfem.h>

#ifdef CALFEM_WITH_MESH

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace calfem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const double pi = 3.14159265358979323846;

// Rectangle w x h with curve markers 10 (bottom), 20 (right), 30 (top), 40 (left)
Geometry rectangle(double w, double h, int nx = -1, int ny = -1)
{
    Geometry g;
    g.point({0.0, 0.0});
    g.point({w, 0.0});
    g.point({w, h});
    g.point({0.0, h});
    g.spline({0, 1}, -1, 10, nx);
    g.spline({1, 2}, -1, 20, ny);
    g.spline({2, 3}, -1, 30, nx);
    g.spline({3, 0}, -1, 40, ny);
    return g;
}

// Rectangle with a circular hole (marker 50) of radius r at (cx, cy)
Geometry plateWithHole(double w, double h, double cx, double cy, double r, double size, double holeSize)
{
    Geometry g;
    g.point({0.0, 0.0}, -1, 0, size);
    g.point({w, 0.0}, -1, 0, size);
    g.point({w, h}, -1, 0, size);
    g.point({0.0, h}, -1, 0, size);
    g.point({cx, cy}, -1, 0, holeSize);           // 4: centre
    g.point({cx + r, cy}, -1, 0, holeSize);       // 5
    g.point({cx, cy + r}, -1, 0, holeSize);       // 6
    g.point({cx - r, cy}, -1, 0, holeSize);       // 7
    g.point({cx, cy - r}, -1, 0, holeSize);       // 8

    g.spline({0, 1}, -1, 10);
    g.spline({1, 2}, -1, 20);
    g.spline({2, 3}, -1, 30);
    g.spline({3, 0}, -1, 40);
    g.circle({5, 4, 6}, -1, 50);
    g.circle({6, 4, 7}, -1, 50);
    g.circle({7, 4, 8}, -1, 50);
    g.circle({8, 4, 5}, -1, 50);

    g.surface({0, 1, 2, 3}, {{4, 5, 6, 7}}, -1, 100);
    return g;
}

double elementArea(const Matrix &coords, const IntMatrix &topo, Eigen::Index e)
{
    double a = 0.0;
    auto n = topo.cols();
    for (Eigen::Index k = 0; k < n; k++)
    {
        auto i = topo(e, k) - 1, j = topo(e, (k + 1) % n) - 1;
        a += coords(i, 0) * coords(j, 1) - coords(j, 0) * coords(i, 1);
    }
    return 0.5 * a;
}

double minAngleDegrees(const Matrix &coords, const IntMatrix &topo, Eigen::Index e)
{
    double minAngle = 180.0;
    auto n = topo.cols();
    for (Eigen::Index k = 0; k < n; k++)
    {
        Eigen::Vector2d p = coords.row(topo(e, k) - 1).transpose();
        Eigen::Vector2d a = coords.row(topo(e, (k + 1) % n) - 1).transpose() - p;
        Eigen::Vector2d b = coords.row(topo(e, (k + n - 1) % n) - 1).transpose() - p;
        minAngle = std::min(minAngle, std::acos(a.dot(b) / (a.norm() * b.norm())) * 180.0 / pi);
    }
    return minAngle;
}

std::vector<int> bdofNodes(const Mesh &m, int marker)
{
    std::vector<int> nodes;
    for (int dof : m.bdofs.at(marker))
        if ((dof - 1) % m.dofs.cols() == 0)
            nodes.push_back((dof - 1) / static_cast<int>(m.dofs.cols()) + 1);
    return nodes;
}

} // namespace

TEST_CASE("Geometry: IDs, markers and queries", "[geometry]")
{
    Geometry g;
    g.point({0.0, 0.0});
    g.point({1.0, 0.0}, -1, 7, 0.5);
    g.point({1.0, 1.0}, 10); // explicit ID: later IDs use the smallest free one
    g.point({0.0, 1.0});

    CHECK(g.points().count(2) == 1); // smallest free ID after 0, 1, 10
    CHECK(g.points().at(1).marker == 7);
    CHECK(g.points().at(1).el_size == 0.5);
    CHECK_FALSE(g.is3D());

    g.spline({0, 1});
    g.circle({1, 0, 2});
    CHECK(g.curves().size() == 2);
    CHECK(g.pointsOnCurves({0, 1}) == std::vector<int>{0, 1, 2});

    g.setCurveMarker(0, 5);
    CHECK(g.curves().at(0).marker == 5);

    auto [xmin, xmax, ymin, ymax] = g.bounding_box_2d();
    CHECK(xmax == 1.0);
    CHECK(ymax == 1.0);

    CHECK_THROWS_AS(g.circle({0, 1}), std::invalid_argument);
    CHECK_THROWS_AS(g.spline({0, 1}, -1, 0, 4, "linear", 1.0), std::invalid_argument);
    CHECK_THROWS_AS(g.setPointMarker(99, 1), std::invalid_argument);

    g.point({0.0, 0.0, 1.0});
    CHECK(g.is3D());
}

TEST_CASE("Built-in mesher: structured rectangle", "[mesh]")
{
    Geometry g = rectangle(3.0, 1.0, 6, 3);
    g.struct_surface({0, 1, 2, 3}, -1, 7);

    SECTION("quadrilaterals")
    {
        MeshGenerator gen(g);
        gen.el_type = 3;
        gen.dofs_per_node = 2;
        gen.backend = MeshBackend::BuiltIn;
        Mesh m = gen.create();

        REQUIRE(m.coords.rows() == 7 * 4);
        REQUIRE(m.topo.rows() == 18);
        REQUIRE(m.edof.cols() == 8);
        CHECK(m.dofs.rows() == 28);
        CHECK(std::all_of(m.element_markers.begin(), m.element_markers.end(), [](int k) { return k == 7; }));

        // Regular grid: all elements 0.5 x 1/3
        for (Eigen::Index e = 0; e < m.topo.rows(); e++)
            CHECK_THAT(elementArea(m.coords, m.topo, e), WithinRel(0.5 / 3.0, 1e-12));

        // Boundary dofs: 7 nodes on the bottom, two dofs each
        CHECK(m.bdofs.at(10).size() == 14);
        for (int node : bdofNodes(m, 10))
            CHECK_THAT(m.coords(node - 1, 1), WithinAbs(0.0, 1e-14));
        for (int node : bdofNodes(m, 20))
            CHECK_THAT(m.coords(node - 1, 0), WithinAbs(3.0, 1e-14));
    }

    SECTION("triangles")
    {
        auto [coords, edof, dofs, bdofs, markers] = create_mesh(g, 2, 1.0, 1);

        CHECK(coords.rows() == 28);
        CHECK(edof.rows() == 36);
        CHECK(edof.cols() == 3);
        CHECK(markers.size() == 36);
        CHECK(bdofs.at(40).size() == 4); // left edge, one dof per node
    }
}

TEST_CASE("Built-in mesher: plate with a hole (triangles)", "[mesh]")
{
    double r = 0.3;
    Geometry g = plateWithHole(2.0, 1.0, 1.0, 0.5, r, 0.1, 0.04);

    MeshGenerator gen(g);
    gen.backend = MeshBackend::BuiltIn;
    gen.el_type = 2;
    Mesh m = gen.create();

    INFO(m.coords.rows() << " nodes, " << m.topo.rows() << " elements");
    REQUIRE(m.topo.rows() > 100);

    // Hole nodes on the circle
    auto holeNodes = bdofNodes(m, 50);
    REQUIRE(holeNodes.size() > 20);
    for (int node : holeNodes)
    {
        double dx = m.coords(node - 1, 0) - 1.0, dy = m.coords(node - 1, 1) - 0.5;
        CHECK_THAT(std::sqrt(dx * dx + dy * dy), WithinRel(r, 1e-12));
    }

    // Elements cover the plate minus the polygonal hole exactly
    double area = 0.0, minArea = 1e300, worstAngle = 180.0;
    for (Eigen::Index e = 0; e < m.topo.rows(); e++)
    {
        double a = elementArea(m.coords, m.topo, e);
        area += a;
        minArea = std::min(minArea, a);
        worstAngle = std::min(worstAngle, minAngleDegrees(m.coords, m.topo, e));
    }

    // Polygon area of the hole from its boundary nodes, in angular order
    std::vector<std::pair<double, int>> byAngle;
    for (int node : holeNodes)
        byAngle.emplace_back(std::atan2(m.coords(node - 1, 1) - 0.5, m.coords(node - 1, 0) - 1.0), node);
    std::sort(byAngle.begin(), byAngle.end());
    double holeArea = 0.0;
    for (size_t i = 0; i < byAngle.size(); i++)
    {
        auto a = byAngle[i].second - 1, b = byAngle[(i + 1) % byAngle.size()].second - 1;
        holeArea += m.coords(a, 0) * m.coords(b, 1) - m.coords(b, 0) * m.coords(a, 1);
    }
    holeArea *= 0.5;

    INFO("min angle " << worstAngle << " degrees");
    CHECK(minArea > 0.0);
    CHECK_THAT(area, WithinRel(2.0 - holeArea, 1e-12));
    CHECK(worstAngle > 20.0);

    // Graded: smaller elements near the hole
    CHECK(holeNodes.size() > 35);
    CHECK(std::all_of(m.element_markers.begin(), m.element_markers.end(), [](int k) { return k == 100; }));
}

TEST_CASE("Built-in mesher: unstructured quadrilaterals", "[mesh]")
{
    Geometry g = plateWithHole(2.0, 1.0, 1.0, 0.5, 0.3, 0.1, 0.05);

    MeshGenerator gen(g);
    gen.backend = MeshBackend::BuiltIn;
    gen.el_type = 3;
    gen.dofs_per_node = 2;
    Mesh m = gen.create();

    REQUIRE(m.topo.cols() == 4);

    double area = 0.0, worstAngle = 180.0;
    for (Eigen::Index e = 0; e < m.topo.rows(); e++)
    {
        double a = elementArea(m.coords, m.topo, e);
        CHECK(a > 0.0);
        area += a;
        worstAngle = std::min(worstAngle, minAngleDegrees(m.coords, m.topo, e));
    }

    INFO("min angle " << worstAngle << " degrees");
    CHECK(area > 2.0 - 0.3 * 0.3 * pi - 0.01);
    CHECK(area < 2.0 - 0.3 * 0.3 * pi + 0.01);
    CHECK(worstAngle > 15.0);

    // Odd el_on_curve is not possible for quadrilaterals on unstructured surfaces
    Geometry odd = rectangle(1.0, 1.0, 5, 4);
    odd.surface({0, 1, 2, 3});
    MeshGenerator oddGen(odd);
    oddGen.backend = MeshBackend::BuiltIn;
    oddGen.el_type = 3;
    CHECK_THROWS_AS(oddGen.create(), std::invalid_argument);
}

TEST_CASE("Built-in mesher: curve types and distributions", "[mesh]")
{
    Geometry g;
    g.point({0.0, 0.0});  // 0 centre
    g.point({2.0, 0.0});  // 1 ellipse start (on the major axis)
    g.point({0.0, 1.0});  // 2 ellipse end
    g.point({-2.0, 0.0}); // 3
    g.point({0.0, -1.0}); // 4

    // Ellipse with a = 2, b = 1 from four quarter arcs
    g.ellipse({1, 0, 1, 2}, -1, 1, 8);
    g.ellipse({2, 0, 1, 3}, -1, 1, 8);
    g.ellipse({3, 0, 1, 4}, -1, 1, 8);
    g.ellipse({4, 0, 1, 1}, -1, 1, 8);
    g.surface({0, 1, 2, 3});

    MeshGenerator gen(g);
    gen.backend = MeshBackend::BuiltIn;
    Mesh m = gen.create();

    auto nodes = bdofNodes(m, 1);
    CHECK(nodes.size() == 32);
    for (int node : nodes)
    {
        double x = m.coords(node - 1, 0), y = m.coords(node - 1, 1);
        CHECK_THAT(x * x / 4.0 + y * y, WithinAbs(1.0, 1e-12));
    }

    // Progression: element lengths grow by the factor
    Geometry p = rectangle(1.0, 1.0);
    p.spline({0, 2}, 10, 5, 4, "progression", 2.0);
    p.surface({0, 1, 2, 3});
    MeshGenerator pgen(p);
    pgen.backend = MeshBackend::BuiltIn;
    pgen.el_size_factor = 0.25;
    Mesh pm = pgen.create();

    std::vector<double> d;
    for (int node : bdofNodes(pm, 5))
        d.push_back(pm.coords(node - 1, 0));
    std::sort(d.begin(), d.end());
    REQUIRE(d.size() == 5);
    CHECK_THAT((d[2] - d[1]) / (d[1] - d[0]), WithinRel(2.0, 1e-9));
    CHECK_THAT((d[4] - d[3]) / (d[3] - d[2]), WithinRel(2.0, 1e-9));
}

TEST_CASE("Built-in mesher: element size factor and shared curves", "[mesh]")
{
    // Two surfaces sharing a curve, one structured and one unstructured
    Geometry g;
    g.point({0.0, 0.0});
    g.point({1.0, 0.0});
    g.point({1.0, 1.0});
    g.point({0.0, 1.0});
    g.point({2.0, 0.0});
    g.point({2.0, 1.0});
    g.spline({0, 1}, -1, 0, 4);
    g.spline({1, 2}, -1, 1, 4); // shared
    g.spline({2, 3}, -1, 0, 4);
    g.spline({3, 0}, -1, 0, 4);
    g.spline({1, 4}, -1, 0);
    g.spline({4, 5}, -1, 0);
    g.spline({5, 2}, -1, 0);
    g.struct_surface({0, 1, 2, 3}, -1, 1);
    g.surface({4, 5, 6, 1}, {}, -1, 2);

    MeshGenerator gen(g);
    gen.backend = MeshBackend::BuiltIn;
    gen.el_size_factor = 0.2;
    Mesh m = gen.create();

    // The shared curve has 5 nodes, used by elements of both surfaces
    auto shared = bdofNodes(m, 1);
    REQUIRE(shared.size() == 5);

    std::set<int> usedBy1, usedBy2;
    for (Eigen::Index e = 0; e < m.topo.rows(); e++)
        for (Eigen::Index k = 0; k < m.topo.cols(); k++)
            (m.element_markers[e] == 1 ? usedBy1 : usedBy2).insert(m.topo(e, k));

    for (int node : shared)
    {
        CHECK(usedBy1.count(node) == 1);
        CHECK(usedBy2.count(node) == 1);
    }

    double area = 0.0;
    for (Eigen::Index e = 0; e < m.topo.rows(); e++)
        area += elementArea(m.coords, m.topo, e);
    CHECK_THAT(area, WithinRel(2.0, 1e-12));
}

TEST_CASE("Built-in mesher: unsupported input", "[mesh]")
{
    Geometry g3;
    g3.point({0.0, 0.0, 0.0});
    g3.point({1.0, 0.0, 0.0});
    g3.point({1.0, 1.0, 0.0});
    g3.spline({0, 1});
    g3.spline({1, 2});
    g3.spline({2, 0});
    g3.surface({0, 1, 2});

    MeshGenerator gen(g3);
    gen.backend = MeshBackend::BuiltIn;
    CHECK_THROWS_AS(gen.create(), std::invalid_argument);

    Geometry g = rectangle(1.0, 1.0);
    g.surface({0, 1, 2, 3});
    MeshGenerator gen2(g);
    gen2.backend = MeshBackend::BuiltIn;
    gen2.el_type = 9;
    CHECK_THROWS_AS(gen2.create(), std::invalid_argument);

    Geometry open = rectangle(1.0, 1.0);
    open.surface({0, 1, 2});
    MeshGenerator gen3(open);
    gen3.backend = MeshBackend::BuiltIn;
    CHECK_THROWS_AS(gen3.create(), std::invalid_argument);
}

TEST_CASE("apply_traction_linear_element: total force", "[mesh][boundary]")
{
    Geometry g = rectangle(2.0, 1.0, 4, 2);
    g.struct_surface({0, 1, 2, 3});

    MeshGenerator gen(g);
    gen.backend = MeshBackend::BuiltIn;
    gen.el_type = 3;
    gen.dofs_per_node = 2;
    Mesh m = gen.create();

    ColVec F = ColVec::Zero(m.dofs.size());
    RowVec q(2);
    q << 0.0, -3.0;
    apply_traction_linear_element(m.boundary_elements, m.coords, m.dofs, F, 30, q); // top edge, length 2

    double Fy = 0.0;
    for (Eigen::Index i = 0; i < m.dofs.rows(); i++)
        Fy += F(m.dofs(i, 1) - 1);
    CHECK_THAT(Fy, WithinRel(-6.0, 1e-12));
}

#endif
