//
// CALFEM C++ example: 3D truss with result plots
//
// The truss of truss3d (three bars from fixed supports to a loaded top node) drawn with
// calfem::vis in an orthographic projection: the undeformed and deformed truss with the supports,
// the load and the normal force of each bar.
//
// Built with Qt (CALFEM_WITH_QT) the figure is shown in a window; otherwise, or with the
// argument --svg, it is written to truss3d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include "example_vis.h"

#include <cstdio>

using namespace calfem;
using namespace example;

int main(int argc, char *argv[])
{
    // ---- Model (as truss3d) ----

    Matrix coords(4, 3);
    coords << 0.0, 0.0, 0.0, //
        2.0, 0.0, 0.0,       //
        0.0, 2.0, 0.0,       //
        0.5, 0.5, 2.0;

    IntMatrix dofs(4, 3);
    dofs << 1, 2, 3, //
        4, 5, 6,     //
        7, 8, 9,     //
        10, 11, 12;

    IntMatrix topo(3, 2);
    topo << 1, 4, //
        2, 4,     //
        3, 4;

    RowVec ep(2);
    ep << 210e9, 1e-4;

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    TripletList triplets;
    for (Eigen::Index i = 0; i < topo.rows(); i++)
        spassem(edof.row(i), triplets, bar3e(ex.row(i), ey.row(i), ez.row(i), ep));

    SpMatrix K(12, 12);
    K.setFromTriplets(triplets.begin(), triplets.end());

    ColVec f = ColVec::Zero(12);
    double P = -10e3;
    f(dofs(3, 2) - 1) += P;

    IntColVec bcDofs(9);
    bcDofs << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    auto [a, r] = spsolveq(K, f, bcDofs);

    Matrix ed = extract_eldisp(edof, a);

    // ---- Projection to the drawing plane ----
    //
    // The projection is linear, so the projected displacements of the bar ends give
    // the deformed bars with eldisp2 (bar elements: ed = [u1, v1, u2, v2]).

    auto nel = topo.rows();
    Matrix px(nel, 2), py(nel, 2), ped(nel, 4);
    for (Eigen::Index i = 0; i < nel; i++)
        for (int n = 0; n < 2; n++)
        {
            vis::Point p = project(ex(i, n), ey(i, n), ez(i, n));
            vis::Point u = project(ed(i, 3 * n), ed(i, 3 * n + 1), ed(i, 3 * n + 2));
            px(i, n) = p.x;
            py(i, n) = p.y;
            ped(i, 2 * n) = u.x;
            ped(i, 2 * n + 1) = u.y;
        }

    // ---- Plot ----

    vis::Figure &fig = vis::figure(1);
    axes3d(fig, -0.5, 2.5, 0.0, 0.5);
    vis::eldraw2(fig, px, py, {1, 1, 1}, IntRowVec::LinSpaced(nel, 1, int(nel)));
    double sfac = vis::eldisp2(fig, px, py, ped, {2, 4, 2});

    // Supports and load
    vis::Markers supports;
    for (int n = 0; n < 3; n++)
        supports.points.push_back(project(coords(n, 0), coords(n, 1), coords(n, 2)));
    supports.marker = vis::Marker::Circle;
    supports.color = colorBlack;
    supports.size = 10.0;
    supports.filled = true;
    fig.add(supports);

    vis::Point top = project(coords(3, 0), coords(3, 1), coords(3, 2));
    arrow(fig, top.x, top.y + 0.6, 0.0, -0.55, colorGreen);
    vis::text(fig, "P = 10 kN", {top.x + 0.08, top.y + 0.35}, colorGreen);

    // Normal forces at the bar midpoints
    for (Eigen::Index i = 0; i < nel; i++)
    {
        ColVec es = bar3s(ex.row(i), ey.row(i), ez.row(i), ep, ed.row(i));
        char label[32];
        std::snprintf(label, sizeof(label), "N = %.2f kN", es(0) / 1e3);
        vis::text(fig, label, {0.5 * (px(i, 0) + px(i, 1)) + 0.05, 0.5 * (py(i, 0) + py(i, 1))}, colorBlue, 9.0);
    }

    char title[64];
    std::snprintf(title, sizeof(title), "Truss, displacements magnified %.0f times", sfac);
    vis::title(fig, title);

    return finish(argc, argv, "truss3d", {"displacements"});
}
