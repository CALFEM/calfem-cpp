//
// CALFEM C++ example: 2D portal frame with result plots
//
// The frame of frame2d (two clamped columns, a beam with a distributed load and a horizontal
// point load) drawn with calfem::vis as in calfem-python's exs_beam2: the undeformed and
// deformed frame, and normal force, shear force and moment diagrams.
//
// Built with Qt (CALFEM_WITH_QT) the figures are shown in windows; otherwise, or with the
// argument --svg, they are written to frame2d_*.svg.
//
// SPDX-License-Identifier: MIT
//

#include <calfem/calfem.h>

#ifdef CALFEM_EXAMPLE_QT
#include <calfem/vis_qt.h>
#endif

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace calfem;

int main(int argc, char *argv[])
{
    // ---- Model (as frame2d) ----

    Matrix coords(4, 2);
    coords << 0.0, 0.0, //
        0.0, 4.0,       //
        6.0, 4.0,       //
        6.0, 0.0;

    IntMatrix dofs(4, 3);
    dofs << 1, 2, 3, //
        4, 5, 6,     //
        7, 8, 9,     //
        10, 11, 12;

    IntMatrix topo(3, 2);
    topo << 1, 2, //
        2, 3,     //
        3, 4;

    RowVec ep(3); // E [Pa], A [m^2], I [m^4]
    ep << 210e9, 45.3e-4, 2510e-8;

    Matrix eq = Matrix::Zero(3, 2);
    eq(1, 1) = -10e3; // 10 kN/m on the beam

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    Matrix K = Matrix::Zero(12, 12);
    ColVec f = ColVec::Zero(12);
    for (Eigen::Index i = 0; i < topo.rows(); i++)
    {
        auto [Ke, fe] = beam2e(ex.row(i), ey.row(i), ep, eq.row(i));
        assem(edof.row(i), K, Ke, f, fe);
    }
    f(dofs(1, 0) - 1) += 2e3;

    IntColVec bcDofs(6);
    bcDofs << 1, 2, 3, 10, 11, 12;
    auto [a, r] = solveq(K, f, bcDofs);

    Matrix ed = extract_eldisp(edof, a);

    // Section forces and local displacements in 21 points per element
    const int nep = 21;
    std::vector<Matrix> es(3), edi(3), eci(3);
    for (Eigen::Index i = 0; i < topo.rows(); i++)
        std::tie(es[i], edi[i], eci[i]) = beam2s(ex.row(i), ey.row(i), ep, ed.row(i), eq.row(i), nep);

    // ---- Plots ----

    std::array<int, 3> plotpar = {2, 1, 0};
    const std::array<double, 4> limits = {-1.5, 7.5, -0.5, 5.5};

    // Undeformed and deformed frame, one scale factor for all elements
    vis::figure(1);
    vis::eldraw2(ex, ey, {1, 2, 1}, IntRowVec::LinSpaced(3, 1, 3));
    double sfac = vis::eldisp2(ex, ey, ed, plotpar);
    vis::axis(limits);
    vis::scalgraph2(sfac, {5e-3, 0.5, 0.0}, 1);
    vis::title("Displacements");

    // Section force diagrams with a common scale factor per quantity
    const char *titles[] = {"Normal force", "Shear force", "Moment"};
    const double scaleBars[] = {3e4, 3e4, 3e4};
    for (int q = 0; q < 3; q++)
    {
        double esMax = 0.0;
        for (const auto &e : es)
            esMax = std::max(esMax, e.col(q).cwiseAbs().maxCoeff());
        double scale = 0.2 * 6.0 / esMax; // 0.2 times the beam length

        vis::figure(q + 2);
        for (Eigen::Index i = 0; i < topo.rows(); i++)
            vis::secforce2(ex.row(i), ey.row(i), es[i].col(q), {2, 1}, scale);
        vis::axis(limits);
        vis::scalgraph2(scale, {scaleBars[q], 0.5, 0.0}, 1);
        vis::title(titles[q]);
    }

    bool svg = argc > 1 && std::string(argv[1]) == "--svg";

#ifdef CALFEM_EXAMPLE_QT
    if (!svg)
    {
        vis::show();
        return 0;
    }
#endif

    const char *names[] = {"displacements", "normal_force", "shear_force", "moment"};
    for (int n = 1; n <= 4; n++)
    {
        std::string filename = std::string("frame2d_") + names[n - 1] + ".svg";
        vis::save_figure(filename, vis::figure(n));
        std::cout << "Figure " << n << " written to " << filename << "\n";
    }
    (void)svg;

    return 0;
}
