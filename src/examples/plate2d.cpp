//
// CALFEM C++ example: simply supported square plate
//
// A 2 m x 2 m steel plate (thickness 20 mm) with a uniform load is simply supported along
// all edges. It is meshed with rectangular Kirchhoff plate elements (platre). The edge dofs
// are collected in a BoundaryDofs map and fixed with apply_bc_3d (three dofs per node). The
// centre deflection is compared with Navier's solution w = 0.00406 q a^4 / D.
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

    double a = 2.0, t = 0.02;
    double E = 210e9, nu = 0.3;
    double q = -10e3; // uniform load [N/m^2]

    int n = 10; // elements per side

    // ---- Mesh, three dofs (w, theta_x, theta_y) per node ----

    int nnodes = (n + 1) * (n + 1);
    auto node = [n](int i, int j) { return j * (n + 1) + i + 1; }; // 1-based

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

    // ---- Assembly ----

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

    // ---- Boundary conditions: w = 0 on all edges (marker 1) ----

    BoundaryDofs bdofs;
    for (int k = 0; k <= n; k++)
        for (int nd : {node(k, 0), node(k, n), node(0, k), node(n, k)})
            for (int d = 0; d < 3; d++)
                bdofs[1].push_back(dofs(nd - 1, d));

    auto [bc, bcVal] = apply_bc_3d(bdofs, IntColVec(), ColVec(), 1, 0.0, 1); // dimension 1: w only

    auto [w, r] = spsolveq(K, f, bc, bcVal);

    // ---- Centre deflection ----

    double Dplate = E * t * t * t / (12.0 * (1.0 - nu * nu));
    double navier = 0.00406235 * q * a * a * a * a / Dplate;
    double centre = w(dofs(node(n / 2, n / 2) - 1, 0) - 1);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nCentre deflection [m]:\n";
    std::cout << "  FE (platre, " << n << " x " << n << " elements): " << centre << "\n";
    std::cout << "  Navier solution:                  " << navier << "\n";
    std::cout << "  Ratio FE / Navier:                " << std::fixed << centre / navier << "\n";

    // ---- Equilibrium: support reactions balance the load ----

    double reactions = 0.0;
    for (Eigen::Index i = 0; i < bc.size(); i++)
        reactions += r(bc(i) - 1);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nSum of support reactions: " << reactions << " N (total load " << q * a * a << " N)\n";
    std::cout << "Prescribed dofs:          " << bc.size() << "\n";

    return 0;
}
