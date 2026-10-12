//
// CALFEM C++ example: simply supported square plate with result plots
//
// The plate of plate2d (platre elements, uniform load, simply supported edges) drawn with
// calfem::vis: the deflected plate as a wireframe (orthographic projection), and the
// deflection along the centre line compared with Navier's series solution.
//
// Built with Qt (CALFEM_WITH_QT) the figures are shown in windows; otherwise, or with the
// argument --svg, they are written to plate2d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include "example_vis.h"

#include <algorithm>
#include <cstdio>

using namespace calfem;
using namespace example;

int main(int argc, char *argv[])
{
    const double pi = 3.14159265358979323846;

    // ---- Model (as plate2d) ----

    double a = 2.0, t = 0.02;
    double E = 210e9, nu = 0.3;
    double q = -10e3;
    int n = 10;

    int nnodes = (n + 1) * (n + 1);
    auto node = [n](int i, int j) { return j * (n + 1) + i + 1; };

    Matrix coords(nnodes, 2);
    for (int j = 0; j <= n; j++)
        for (int i = 0; i <= n; i++)
            coords.row(node(i, j) - 1) << a * i / n, a * j / n;

    IntMatrix dofs = create_dofs(nnodes, 3);

    IntMatrix topo(n * n, 4);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            topo.row(j * n + i) << node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1);

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    RowVec ep(1);
    ep << t;
    Matrix D = hooke(1, E, nu);

    int ndof = 3 * nnodes;
    TripletList triplets;
    ColVec f = ColVec::Zero(ndof);
    for (Eigen::Index e = 0; e < topo.rows(); e++)
    {
        auto [Ke, fe] = platre(ex.row(e), ey.row(e), ep, D, q);
        spassem(edof.row(e), triplets, Ke, f, fe);
    }

    SpMatrix K(ndof, ndof);
    K.setFromTriplets(triplets.begin(), triplets.end());

    BoundaryDofs bdofs;
    for (int k = 0; k <= n; k++)
        for (int nd : {node(k, 0), node(k, n), node(0, k), node(n, k)})
            for (int d = 0; d < 3; d++)
                bdofs[1].push_back(dofs(nd - 1, d));

    auto [bc, bcVal] = apply_bc_3d(bdofs, IntColVec(), ColVec(), 1, 0.0, 1);
    auto [w, r] = spsolveq(K, f, bc, bcVal);

    // ---- Figure 1: deflected plate in an orthographic projection ----
    //
    // The edges stay at w = 0, so the boundary of the deflected mesh shows the undeformed plate.

    double wmax = 0.0;
    for (int nd = 0; nd < nnodes; nd++)
        wmax = std::max(wmax, std::abs(w(dofs(nd, 0) - 1)));
    double sfac = 0.4 / wmax; // largest deflection drawn as 0.4 m

    auto nel = topo.rows();
    Matrix px(nel, 4), py(nel, 4);
    for (Eigen::Index e = 0; e < nel; e++)
        for (int k = 0; k < 4; k++)
        {
            int nd = topo(e, k) - 1;
            double x = coords(nd, 0), y = coords(nd, 1);
            vis::Point p = project(x, y, sfac * w(dofs(nd, 0) - 1));
            px(e, k) = p.x;
            py(e, k) = p.y;
        }

    vis::Figure &surface = vis::figure(1);
    axes3d(surface, -0.5, -0.5, 0.0, 0.5);
    vis::eldraw2(surface, px, py, {1, 2, 0});

    char title[80];
    std::snprintf(title, sizeof(title), "Deflected plate, max %.2f mm (magnified %.0f times)", wmax * 1e3, sfac);
    vis::title(surface, title);

    // ---- Figure 2: deflection along y = a/2, FE (markers) and Navier (line) ----

    double Dplate = E * t * t * t / (12.0 * (1.0 - nu * nu));
    auto navier = [&](double x, double y) {
        double sum = 0.0;
        for (int mi = 1; mi <= 39; mi += 2)
            for (int ni = 1; ni <= 39; ni += 2)
                sum += std::sin(mi * pi * x / a) * std::sin(ni * pi * y / a) /
                       (mi * ni * std::pow(mi * mi + ni * ni, 2));
        return 16.0 * q * std::pow(a, 4) / (std::pow(pi, 6) * Dplate) * sum;
    };

    std::vector<double> xn, wn, xs, ws;
    for (int i = 0; i <= n; i++)
    {
        int nd = node(i, n / 2) - 1;
        xn.push_back(coords(nd, 0));
        wn.push_back(w(dofs(nd, 0) - 1) * 1e3);
    }
    for (int i = 0; i <= 100; i++)
    {
        xs.push_back(a * i / 100.0);
        ws.push_back(navier(xs.back(), a / 2.0) * 1e3);
    }

    vis::Figure &centre = vis::figure(2);
    graph(centre);
    plot(centre, xs, ws, colorBlue);
    plotMarkers(centre, xn, wn, colorRed);
    vis::title(centre, "Deflection [mm] along y = a/2: FE (markers), Navier (line)");

    return finish(argc, argv, "plate2d", {"surface", "centre_line"});
}
