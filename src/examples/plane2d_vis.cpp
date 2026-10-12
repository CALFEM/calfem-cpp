//
// CALFEM C++ example: plane stress cantilever with result plots
//
// The cantilever of plane2d (plani4e elements, clamped left edge, tip load) drawn with
// calfem::vis: the undeformed and deformed mesh, the deflection of the centre line
// compared with beam theory including shear deformation, and the effective stress.
//
// Built with Qt (CALFEM_WITH_QT) the figures are shown in windows; otherwise, or with the
// argument --svg, they are written to plane2d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include "example_vis.h"

#include <cstdio>

using namespace calfem;
using namespace example;

int main(int argc, char *argv[])
{
    // ---- Model (as plane2d) ----

    double L = 2.0, h = 0.2, t = 0.01;
    double E = 210e9, nu = 0.3;
    double P = -10e3;
    int nx = 40, ny = 8;

    int nnodes = (nx + 1) * (ny + 1);
    auto node = [nx](int i, int j) { return j * (nx + 1) + i + 1; };

    Matrix coords(nnodes, 2);
    IntMatrix dofs(nnodes, 2);
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
        {
            int n = node(i, j) - 1;
            coords.row(n) << L * i / nx, h * j / ny - h / 2.0;
            dofs.row(n) << 2 * n + 1, 2 * n + 2;
        }

    IntMatrix topo(nx * ny, 4);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            topo.row(j * nx + i) << node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1);

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    int ptype = 1;
    RowVec ep(3);
    ep << ptype, t, 2;
    Matrix D = hooke(ptype, E, nu);
    int ndof = 2 * nnodes;

    TripletList triplets;
    for (Eigen::Index e = 0; e < topo.rows(); e++)
        spassem(edof.row(e), triplets, plani4e(ex.row(e), ey.row(e), ep, D));

    SpMatrix K(ndof, ndof);
    K.setFromTriplets(triplets.begin(), triplets.end());

    ColVec f = ColVec::Zero(ndof);
    for (int j = 0; j <= ny; j++)
        f(dofs(node(nx, j) - 1, 1) - 1) = P / ny * ((j == 0 || j == ny) ? 0.5 : 1.0);

    IntColVec bcDofs(2 * (ny + 1));
    for (int j = 0; j <= ny; j++)
    {
        bcDofs(2 * j) = dofs(node(0, j) - 1, 0);
        bcDofs(2 * j + 1) = dofs(node(0, j) - 1, 1);
    }

    auto [a, r] = spsolveq(K, f, bcDofs);
    Matrix ed = extract_eldisp(edof, a);

    // ---- Figure 1: undeformed and deformed mesh ----

    vis::Figure &mesh = vis::figure(1);
    vis::eldraw2(mesh, ex, ey, {3, 1, 0});
    double sfac = vis::eldisp2(mesh, ex, ey, ed, {1, 2, 0}, 10.0);
    vis::axis(mesh, {-0.1, 2.1, -0.4, 0.2});
    vis::scalgraph2(mesh, sfac, {0.01, 0.0, -0.32}, 1);

    char title[64];
    std::snprintf(title, sizeof(title), "Deformed cantilever (displacements x %.0f)", sfac);
    vis::title(mesh, title);

    // ---- Figure 2: centre line deflection, FE (markers) and beam theory (line) ----

    double I = t * h * h * h / 12.0, A = t * h;
    double G = E / (2.0 * (1.0 + nu));
    auto beam = [&](double x) { return P / (E * I) * (L * x * x / 2.0 - x * x * x / 6.0) + P * x / (5.0 / 6.0 * G * A); };

    std::vector<double> xn, wn, xb, wb;
    for (int i = 0; i <= nx; i++)
    {
        int n = node(i, ny / 2) - 1;
        xn.push_back(coords(n, 0));
        wn.push_back(a(dofs(n, 1) - 1) * 1e3);
    }
    for (int i = 0; i <= 100; i++)
    {
        xb.push_back(L * i / 100.0);
        wb.push_back(beam(xb.back()) * 1e3);
    }

    vis::Figure &deflection = vis::figure(2);
    graph(deflection);
    plot(deflection, xb, wb, colorBlue);
    plotMarkers(deflection, xn, wn, colorRed);
    vis::title(deflection, "Deflection [mm] along x [m]: FE (markers), beam theory (line)");

    // ---- Figure 3: effective stress at the element centres on the deformed mesh ----

    RowVec ep1 = ep;
    ep1(2) = 1; // one Gauss point: the element centre

    ColVec eseff(topo.rows());
    for (Eigen::Index e = 0; e < topo.rows(); e++)
    {
        auto [es, et] = plani4s(ex.row(e), ey.row(e), ep1, D, ed.row(e));
        eseff(e) = effmises(es, ptype)(0) / 1e6;
    }

    vis::Figure &stress = vis::figure(3);
    vis::draw_element_values(stress, eseff, coords, edof, 2, 3, a, false, true, 10.0,
                             "Effective stress [MPa] (displacements x 10)");
    vis::colorbar(stress);
    vis::axis(stress, {-0.1, 2.1, -0.4, 0.2});

    return finish(argc, argv, "plane2d", {"deformed", "deflection", "stress"});
}
