//
// CALFEM C++ example: plate with a hole with result plots
//
// The plate of mesh2d (meshed from a Geometry, pulled at the right edge) drawn with
// calfem::vis: the mesh, the deformed mesh, sigma_xx above the hole compared with
// Kirsch's solution for an infinite plate (the finite width raises the FE stresses somewhat)
// and sigma_xx contours.
//
// Uses the built-in mesher, or Gmsh when calfem is built with CALFEM_WITH_GMSH.
// Arguments (optional): --svg, element type (3 quadrilaterals, default, or 2 triangles).
// Built with Qt (CALFEM_WITH_QT) the figures are shown in windows; otherwise, or with --svg,
// they are written to mesh2d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include "example_vis.h"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>

using namespace calfem;
using namespace example;

int run(int argc, char *argv[])
{
    // ---- Geometry and mesh (as mesh2d) ----

    double w = 0.4, h = 0.2, r = 0.02;
    double cx = w / 2.0, cy = h / 2.0;

    Geometry g;
    g.point({0.0, 0.0}, -1, 0, 0.02);
    g.point({w, 0.0}, -1, 0, 0.02);
    g.point({w, h}, -1, 0, 0.02);
    g.point({0.0, h}, -1, 0, 0.02);

    double hs = 0.003;
    g.point({cx, cy}, -1, 0, hs);
    g.point({cx + r, cy}, -1, 0, hs);
    g.point({cx, cy + r}, -1, 0, hs);
    g.point({cx - r, cy}, -1, 0, hs);
    g.point({cx, cy - r}, -1, 0, hs);

    const int bottom = 10, right = 20, top = 30, left = 40, hole = 50;

    g.spline({0, 1}, -1, bottom);
    g.spline({1, 2}, -1, right);
    g.spline({2, 3}, -1, top);
    g.spline({3, 0}, -1, left);
    g.circle({5, 4, 6}, -1, hole);
    g.circle({6, 4, 7}, -1, hole);
    g.circle({7, 4, 8}, -1, hole);
    g.circle({8, 4, 5}, -1, hole);

    g.surface({0, 1, 2, 3}, {{4, 5, 6, 7}});

    int elType = 3;
    for (int i = 1; i < argc; i++)
        if (std::string(argv[i]) != "--svg")
            elType = std::atoi(argv[i]);

    MeshGenerator mesher(g);
    mesher.el_type = elType;
    mesher.dofs_per_node = 2;
    Mesh m = mesher.create();
    bool quads = mesher.el_type == 3;

    // ---- Plane stress analysis (as mesh2d) ----

    double E = 210e9, nu = 0.3, t = 0.01;
    double sigma0 = 100e6;

    RowVec ep(3);
    ep << 1, t, 2;
    RowVec ept = ep.head(2);
    Matrix D = hooke(1, E, nu);

    auto [ex, ey, ez] = coord_extract(m.edof, m.coords, m.dofs);

    int ndof = static_cast<int>(m.dofs.size());
    TripletList triplets;
    for (Eigen::Index e = 0; e < m.edof.rows(); e++)
    {
        Matrix Ke = quads ? plani4e(ex.row(e), ey.row(e), ep, D) : plante(ex.row(e), ey.row(e), ept, D);
        spassem(m.edof.row(e), triplets, Ke);
    }

    SpMatrix K(ndof, ndof);
    K.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::Index corner;
    m.coords.rowwise().squaredNorm().minCoeff(&corner);

    auto [bc, bcVal] = apply_bc(m.bdofs, IntColVec(), ColVec(), left, 0.0, 1);
    IntColVec bcAll(bc.size() + 1);
    bcAll << bc, m.dofs(corner, 1);
    ColVec bcValAll = ColVec::Zero(bcAll.size());

    ColVec f = ColVec::Zero(ndof);
    RowVec q(2);
    q << sigma0 * t, 0.0;
    apply_traction_linear_element(m.boundary_elements, m.coords, m.dofs, f, right, q);

    auto [a, rf] = spsolveq(K, f, bcAll, bcValAll);
    Matrix ed = extract_eldisp(m.edof, a);

    // ---- Figure 1: mesh ----

    vis::Figure &mesh = vis::figure(1);
    vis::eldraw2(mesh, ex, ey, {1, 2, 0});
    char title[80];
    std::snprintf(title, sizeof(title), "Mesh: %d nodes, %d %s", int(m.coords.rows()), int(m.edof.rows()),
                  quads ? "quadrilaterals" : "triangles");
    vis::title(mesh, title);

    // ---- Figure 2: deformed mesh ----

    vis::Figure &deformed = vis::figure(2);
    double sfac = vis::eldisp2(deformed, ex, ey, ed, {1, 2, 0}, 200.0);
    std::snprintf(title, sizeof(title), "Deformed plate (displacements x %.0f)", sfac);
    vis::title(deformed, title);

    // ---- Figure 3: sigma_xx above the hole (element centres near x = cx) ----

    RowVec ep1 = ep;
    ep1(2) = 1; // one Gauss point: the element centre

    std::vector<double> yf, sf;
    for (Eigen::Index e = 0; e < m.edof.rows(); e++)
    {
        double xc = ex.row(e).mean(), yc = ey.row(e).mean();
        if (std::abs(xc - cx) > 0.6 * hs || yc < cy)
            continue;
        auto [es, et] = quads ? plani4s(ex.row(e), ey.row(e), ep1, D, ed.row(e))
                              : plants(ex.row(e), ey.row(e), ept, D, ed.row(e));
        yf.push_back((yc - cy) * 1e3);
        sf.push_back(es(0, 0) / sigma0);
    }

    // Kirsch: sigma_xx(0, y) = sigma0 (1 + r^2 / (2 y^2) + 3 r^4 / (2 y^4)) for y >= r
    std::vector<double> yk, sk;
    for (int i = 0; i <= 100; i++)
    {
        double y = r + (cy - r) * i / 100.0;
        yk.push_back(y * 1e3);
        sk.push_back(1.0 + r * r / (2.0 * y * y) + 1.5 * std::pow(r / y, 4));
    }

    vis::Figure &stress = vis::figure(3);
    graph(stress);
    plot(stress, yk, sk, colorBlue);
    plotMarkers(stress, yf, sf, colorRed);
    vis::title(stress, "sigma_xx / sigma_0 above the hole vs distance [mm]: FE (markers), Kirsch (line)");

    // ---- Figure 4: sigma_xx contours (element centre values averaged at the nodes) ----

    ColVec nodal = ColVec::Zero(m.coords.rows()), count = ColVec::Zero(m.coords.rows());
    for (Eigen::Index e = 0; e < m.edof.rows(); e++)
    {
        auto [es, et] = quads ? plani4s(ex.row(e), ey.row(e), ep1, D, ed.row(e))
                              : plants(ex.row(e), ey.row(e), ept, D, ed.row(e));
        for (Eigen::Index k = 0; k < m.edof.cols(); k += 2)
        {
            int node = (m.edof(e, k) - 1) / 2;
            nodal(node) += es(0, 0) / sigma0;
            count(node) += 1.0;
        }
    }
    nodal = nodal.cwiseQuotient(count.cwiseMax(1.0));

    vis::Figure &contours = vis::figure(4);
    vis::draw_nodal_values_contourf(contours, nodal, m.coords, m.edof, 12, "sigma_xx / sigma_0", 2, mesher.el_type);
    vis::colorbar(contours);

    return finish(argc, argv, "mesh2d", {"mesh", "deformed", "stress", "contours"});
}

int main(int argc, char *argv[])
{
    try
    {
        return run(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
