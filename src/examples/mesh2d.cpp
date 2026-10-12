//
// CALFEM C++ example: plate with a hole, meshed from a geometry
//
// A plate with a circular hole is described with Geometry (points, splines, circle arcs,
// a surface with a hole and boundary markers), meshed with create_mesh/MeshGenerator and
// solved in plane stress. The plate is pulled at the right edge; the stress concentration at the
// hole is compared with Kirsch's solution for an infinite plate (factor 3).
//
// Uses the built-in mesher, or Gmsh when calfem is built with CALFEM_WITH_GMSH.
// Arguments (optional): SVG file for a picture of the mesh, element type (3 quadrilaterals,
// default, or 2 triangles):  mesh2d mesh.svg 2
//
// SPDX-License-Identifier: MIT
//

#include <calfem/calfem.h>

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace calfem;

int run(int argc, char *argv[])
{
    // ---- Geometry: 0.4 m x 0.2 m plate, hole radius 0.02 m at the centre ----

    double w = 0.4, h = 0.2, r = 0.02;
    double cx = w / 2.0, cy = h / 2.0;

    Geometry g;
    g.point({0.0, 0.0}, -1, 0, 0.02);
    g.point({w, 0.0}, -1, 0, 0.02);
    g.point({w, h}, -1, 0, 0.02);
    g.point({0.0, h}, -1, 0, 0.02);

    double hs = 0.003; // element size at the hole
    g.point({cx, cy}, -1, 0, hs);     // 4 centre
    g.point({cx + r, cy}, -1, 0, hs); // 5
    g.point({cx, cy + r}, -1, 0, hs); // 6
    g.point({cx - r, cy}, -1, 0, hs); // 7
    g.point({cx, cy - r}, -1, 0, hs); // 8

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

    // ---- Mesh: 4-node quadrilaterals, two dofs per node ----

    MeshGenerator mesher(g);
    mesher.el_type = argc > 2 ? std::atoi(argv[2]) : 3;
    mesher.dofs_per_node = 2;
    Mesh m = mesher.create();

    std::cout << "Mesh (" << (mesher.resolvedBackend() == MeshBackend::Gmsh ? "Gmsh" : "built-in")
              << "): " << m.coords.rows() << " nodes, " << m.topo.rows() << " elements\n";

    if (argc > 1)
    {
        auto [ex, ey, ez] = coord_extract(m.edof, m.coords, m.dofs);
        vis::Figure fig;
        vis::eldraw2(fig, ex, ey, {1, 2, 0});
        vis::title(fig, "Plate with a hole");
        fig.set_size(900.0, 500.0);
        vis::save_figure(argv[1], fig);
        std::cout << "Mesh written to " << argv[1] << "\n";
    }

    // ---- Plane stress analysis ----

    double E = 210e9, nu = 0.3, t = 0.01;
    double sigma0 = 100e6; // applied tension [Pa]

    bool quads = mesher.el_type == 3;

    RowVec ep(3);
    ep << 1, t, 2; // plane stress, thickness, 2 x 2 Gauss points (quadrilaterals)
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

    // Left edge: ux = 0; the bottom-left corner also uy = 0
    Eigen::Index corner;
    m.coords.rowwise().squaredNorm().minCoeff(&corner);

    auto [bc, bcVal] = apply_bc(m.bdofs, IntColVec(), ColVec(), left, 0.0, 1);
    IntColVec bcAll(bc.size() + 1);
    bcAll << bc, m.dofs(corner, 1);
    ColVec bcValAll = ColVec::Zero(bcAll.size());

    // Uniform traction on the right edge
    ColVec f = ColVec::Zero(ndof);
    RowVec q(2);
    q << sigma0 * t, 0.0;
    apply_traction_linear_element(m.boundary_elements, m.coords, m.dofs, f, right, q);

    auto [a, rf] = spsolveq(K, f, bcAll, bcValAll);

    // ---- Stress concentration: largest sigma_xx (Gauss points / constant strain triangles) ----

    Matrix ed = extract_eldisp(m.edof, a);
    double maxSxx = 0.0;
    for (Eigen::Index e = 0; e < m.edof.rows(); e++)
    {
        auto [es, et] = quads ? plani4s(ex.row(e), ey.row(e), ep, D, ed.row(e))
                              : plants(ex.row(e), ey.row(e), ept, D, ed.row(e));
        maxSxx = std::max(maxSxx, es.col(0).maxCoeff());
    }

    // Kirsch's factor 3 for an infinite plate; finite width (d/W = 0.2) gives about 3.14 (Heywood)
    double dW = 2.0 * r / h;
    double heywood = 2.0 + std::pow(1.0 - dW, 3.0);
    double ktNet = heywood; // based on the net section stress
    double ktGross = ktNet / (1.0 - dW);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nStress concentration at the hole (max sigma_xx / applied stress):\n";
    std::cout << "  FE (" << (quads ? "quadrilaterals" : "triangles     ") << "):              " << maxSxx / sigma0
              << "\n";
    std::cout << "  Kirsch, infinite plate:            3.00\n";
    std::cout << "  Heywood, finite width (gross):     " << ktGross << "\n";

    return 0;
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
