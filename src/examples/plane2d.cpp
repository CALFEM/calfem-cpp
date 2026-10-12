//
// CALFEM C++ example: plane stress cantilever
//
// A 2 m x 0.2 m steel cantilever (thickness 0.01 m) is clamped at the left edge and loaded by
// a vertical force at the free end. It is meshed with 4-node isoparametric plane stress
// elements (plani4e). The tip deflection is compared with beam theory including shear
// deformation, and the effective (von Mises) stress is computed with plani4s, effmises and
// stress2nodal.
//
// SPDX-License-Identifier: MIT
//

#include <calfem/calfem.h>

#include <iomanip>
#include <iostream>

using namespace calfem;

int main()
{
    // ---- Problem ----

    double L = 2.0, h = 0.2, t = 0.01;
    double E = 210e9, nu = 0.3;
    double P = -10e3; // total tip load [N]

    int nx = 40, ny = 8;

    // ---- Structured mesh, two dofs (ux, uy) per node ----

    int nnodes = (nx + 1) * (ny + 1);
    auto node = [nx](int i, int j) { return j * (nx + 1) + i + 1; }; // 1-based

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

    // ---- Assembly ----

    int ptype = 1; // plane stress
    RowVec ep(3);
    ep << ptype, t, 2; // ptype, thickness, 2 x 2 Gauss points

    Matrix D = hooke(ptype, E, nu);
    int ndof = 2 * nnodes;

    TripletList triplets;
    for (Eigen::Index e = 0; e < topo.rows(); e++)
        spassem(edof.row(e), triplets, plani4e(ex.row(e), ey.row(e), ep, D));

    SpMatrix K(ndof, ndof);
    K.setFromTriplets(triplets.begin(), triplets.end());

    // ---- Loads and boundary conditions ----

    // Tip load distributed over the nodes of the free edge (half at the corners)
    ColVec f = ColVec::Zero(ndof);
    for (int j = 0; j <= ny; j++)
        f(dofs(node(nx, j) - 1, 1) - 1) = P / ny * ((j == 0 || j == ny) ? 0.5 : 1.0);

    IntColVec bcDofs(2 * (ny + 1)); // clamped left edge
    for (int j = 0; j <= ny; j++)
    {
        bcDofs(2 * j) = dofs(node(0, j) - 1, 0);
        bcDofs(2 * j + 1) = dofs(node(0, j) - 1, 1);
    }

    auto [a, r] = spsolveq(K, f, bcDofs);

    // ---- Tip deflection ----

    double I = t * h * h * h / 12.0, A = t * h;
    double G = E / (2.0 * (1.0 + nu));
    double bending = P * L * L * L / (3.0 * E * I);
    double shear = P * L / (5.0 / 6.0 * G * A);

    double tip = a(dofs(node(nx, ny / 2) - 1, 1) - 1);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nTip deflection [m]:\n";
    std::cout << "  FE (plani4e, " << nx << " x " << ny << " elements): " << tip << "\n";
    std::cout << "  Beam theory (bending + shear):     " << bending + shear << "\n";
    std::cout << "  Ratio FE / beam theory:            " << std::fixed << tip / (bending + shear) << "\n";

    // ---- Stresses ----

    Matrix ed = extract_eldisp(edof, a);

    // Effective stress at the element centres (one Gauss point)
    RowVec ep1 = ep;
    ep1(2) = 1;

    ColVec eseff(topo.rows());
    double maxSxx = 0.0;
    for (Eigen::Index e = 0; e < topo.rows(); e++)
    {
        auto [es, et] = plani4s(ex.row(e), ey.row(e), ep1, D, ed.row(e));
        eseff(e) = effmises(es, ptype)(0);
        maxSxx = std::max(maxSxx, std::abs(es(0, 0)));
    }

    // Nodal values for plotting (each node collects eseff / 4 from its elements)
    Matrix ev = stress2nodal(eseff, edof);

    double sigmaBeam = std::abs(P) * L * (h / 2.0) / I;

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "\nStresses [Pa]:\n";
    std::cout << "  max |sigma_xx| at element centres:  " << maxSxx << "\n";
    std::cout << "  beam theory at the clamped edge:    " << sigmaBeam << " (at the outer fibre)\n";
    std::cout << "  max effective stress (elements):    " << eseff.maxCoeff() << "\n";
    std::cout << "  nodal values from stress2nodal:     " << ev.rows() << " x " << ev.cols() << "\n";

    // ---- Equilibrium ----

    double Ry = 0.0;
    for (Eigen::Index i = 0; i < bcDofs.size(); i += 2)
        Ry += r(bcDofs(i + 1) - 1);

    std::cout << "\nVertical support reaction: " << Ry << " N (applied " << P << " N)\n";

    return 0;
}
