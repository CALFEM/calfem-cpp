//
// CALFEM C++ example: steady-state heat conduction with result plots
//
// The plate of heat2d (uniform heat generation, 0 C on the left edge and 100 C on the right
// edge) drawn with calfem::vis: the mesh with the heat flux at the element centres (flw2i4s),
// and the temperature along the plate compared with the exact solution.
//
// Built with Qt (CALFEM_WITH_QT) the figures are shown in windows; otherwise, or with the
// argument --svg, they are written to heat2d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include "example_vis.h"

#include <algorithm>

using namespace calfem;
using namespace example;

int main(int argc, char *argv[])
{
    // ---- Model (as heat2d) ----

    double Lx = 2.0, Ly = 1.0;
    double k = 1.5, Q = 20.0, T1 = 100.0;
    int nx = 20, ny = 4;

    int nnodes = (nx + 1) * (ny + 1);
    auto node = [nx](int i, int j) { return j * (nx + 1) + i + 1; };

    Matrix coords(nnodes, 2);
    IntMatrix dofs(nnodes, 1);
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
        {
            int n = node(i, j) - 1;
            coords.row(n) << Lx * i / nx, Ly * j / ny;
            dofs(n, 0) = n + 1;
        }

    IntMatrix topo(nx * ny, 4);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            topo.row(j * nx + i) << node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1);

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    RowVec ep(2);
    ep << 1.0, 2;
    Matrix D = k * Matrix::Identity(2, 2);

    TripletList triplets;
    ColVec f = ColVec::Zero(nnodes);
    for (Eigen::Index e = 0; e < topo.rows(); e++)
    {
        auto [Ke, fe] = flw2i4e(ex.row(e), ey.row(e), ep, D, Q);
        spassem(edof.row(e), triplets, Ke, f, fe);
    }

    SpMatrix K(nnodes, nnodes);
    K.setFromTriplets(triplets.begin(), triplets.end());

    IntColVec bcDofs(2 * (ny + 1));
    ColVec bcVals(2 * (ny + 1));
    for (int j = 0; j <= ny; j++)
    {
        bcDofs(2 * j) = node(0, j);
        bcVals(2 * j) = 0.0;
        bcDofs(2 * j + 1) = node(nx, j);
        bcVals(2 * j + 1) = T1;
    }

    auto [a, r] = spsolveq(K, f, bcDofs, bcVals);
    Matrix ed = extract_eldisp(edof, a);

    // ---- Figure 1: mesh and heat flux at the element centres ----

    RowVec ep1 = ep;
    ep1(1) = 1; // one Gauss point: the element centre

    Matrix q(topo.rows(), 2), qc(topo.rows(), 2);
    for (Eigen::Index e = 0; e < topo.rows(); e++)
    {
        auto [es, et, eci] = flw2i4s(ex.row(e), ey.row(e), ep1, D, ed.row(e));
        q.row(e) = es.row(0);
        qc.row(e) = eci.row(0);
    }

    vis::Figure &flux = vis::figure(1);
    vis::eldraw2(flux, ex, ey, {1, 1, 0});

    double qmax = q.rowwise().norm().maxCoeff();
    double scale = 0.09 / qmax; // longest arrow 0.09 m (element width 0.1 m)
    for (Eigen::Index e = 0; e < topo.rows(); e++)
        arrow(flux, qc(e, 0) - 0.5 * scale * q(e, 0), qc(e, 1) - 0.5 * scale * q(e, 1), scale * q(e, 0),
              scale * q(e, 1), colorRed);

    vis::text(flux, "T = 0 C", {0.02, -0.12}, colorBlue);
    vis::text(flux, "T = 100 C", {Lx - 0.35, -0.12}, colorBlue);
    vis::title(flux, "Heat flux at the element centres");

    // ---- Figure 2: temperature along y = 0, FE (markers) and exact (line) ----

    auto exact = [&](double x) { return Q / (2.0 * k) * x * (Lx - x) + T1 * x / Lx; };

    std::vector<double> xn, Tn, xe, Te;
    for (int i = 0; i <= nx; i++)
    {
        xn.push_back(coords(node(i, 0) - 1, 0));
        Tn.push_back(a(node(i, 0) - 1));
    }
    for (int i = 0; i <= 100; i++)
    {
        xe.push_back(Lx * i / 100.0);
        Te.push_back(exact(xe.back()));
    }

    vis::Figure &temperature = vis::figure(2);
    graph(temperature);
    plot(temperature, xe, Te, colorBlue);
    plotMarkers(temperature, xn, Tn, colorRed);
    vis::title(temperature, "Temperature [C] along x [m]: FE (markers), exact (line)");

    return finish(argc, argv, "heat2d", {"flux", "temperature"});
}
