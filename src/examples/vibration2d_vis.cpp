//
// CALFEM C++ example: vibration of a cantilever beam with result plots
//
// The cantilever of vibration2d (10 beam2de elements) drawn with calfem::vis: the first three
// bending mode shapes from eigen, and the tip displacement history from step2 for a short load
// pulse at the free end.
//
// Built with Qt (CALFEM_WITH_QT) the figures are shown in windows; otherwise, or with the
// argument --svg, they are written to vibration2d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include "example_vis.h"

#include <cstdio>

using namespace calfem;
using namespace example;

int main(int argc, char *argv[])
{
    const double pi = 3.14159265358979323846;

    // ---- Model (as vibration2d) ----

    int nel = 10;
    double L = 3.0;
    double E = 210e9, A = 45.3e-4, I = 2510e-8, m = 35.6;

    int nnodes = nel + 1;
    int ndof = 3 * nnodes;

    Matrix coords = Matrix::Zero(nnodes, 2);
    IntMatrix dofs(nnodes, 3);
    for (int i = 0; i < nnodes; i++)
    {
        coords(i, 0) = L * i / nel;
        dofs.row(i) << 3 * i + 1, 3 * i + 2, 3 * i + 3;
    }

    IntMatrix topo(nel, 2);
    for (int e = 0; e < nel; e++)
        topo.row(e) << e + 1, e + 2;

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    RowVec ep(4);
    ep << E, A, I, m;

    Matrix K = Matrix::Zero(ndof, ndof);
    Matrix M = Matrix::Zero(ndof, ndof);
    for (int e = 0; e < nel; e++)
    {
        auto [Ke, Me, Ce] = beam2de(ex.row(e), ey.row(e), ep);
        assem(edof.row(e), K, Ke);
        assem(edof.row(e), M, Me);
    }

    IntColVec bcDofs(3);
    bcDofs << 1, 2, 3;

    auto [lambda, X] = eigen(K, M, bcDofs);

    // ---- Figure 1: bending mode shapes ----

    vis::Figure &modes = vis::figure(1);
    vis::eldraw2(modes, ex, ey, {3, 1, 0});

    const int colors[] = {2, 3, 4};
    const vis::Color labelColors[] = {colorBlue, {0.75, 0.0, 0.75}, colorRed};
    int mode = 0;
    for (Eigen::Index i = 0; (i < lambda.size()) && (mode < 3); i++)
    {
        if (std::abs(X(ndof - 2, i)) < 1e-12) // axial mode
            continue;

        // Mode normalised to a unit tip displacement, drawn with amplitude 0.5 m
        ColVec x = X.col(i) / X(ndof - 2, i);
        Matrix ed = extract_eldisp(edof, x);
        vis::eldisp2(modes, ex, ey, ed, {1, colors[mode], 0}, 0.5);

        char label[48];
        std::snprintf(label, sizeof(label), "mode %d: %.2f Hz", mode + 1, std::sqrt(lambda(i)) / (2.0 * pi));
        vis::text(modes, label, {0.1, 0.7 - 0.12 * mode}, labelColors[mode]);
        mode++;
    }
    vis::axis(modes, {-0.2, 3.4, -0.6, 0.8});
    vis::title(modes, "Bending modes (unit tip displacement)");

    // ---- Figure 2: tip response to a load pulse ----

    double dt = 0.0005, T = 0.2;

    Matrix G(4, 2);
    G << 0.0, 0.0, //
        0.005, 1e3, //
        0.01, 0.0, //
        T, 0.0;

    auto [t, g] = gfunc(G, dt);

    Matrix f = Matrix::Zero(ndof, t.size());
    f.row(ndof - 2) = g.transpose();

    Matrix bc(3, 2);
    bc << 1, 0.0, //
        2, 0.0,   //
        3, 0.0;

    IntColVec tip(1);
    tip << ndof - 1;

    ColVec a0 = ColVec::Zero(ndof);
    RowVec ip(4);
    ip << dt, T, 0.25, 0.5;

    auto [modelhist, dofhist] = step2(K, Matrix(), M, f, a0, a0, bc, ip, RowVec(), tip);

    std::vector<double> time, uy;
    for (Eigen::Index k = 0; k < dofhist.a.cols(); k++)
    {
        time.push_back(k * dt);
        uy.push_back(dofhist.a(0, k) * 1e3);
    }

    vis::Figure &response = vis::figure(2);
    graph(response);
    plot(response, time, uy, colorBlue);
    vis::title(response, "Tip displacement uy [mm] versus time [s]");

    return finish(argc, argv, "vibration2d", {"modes", "response"});
}
