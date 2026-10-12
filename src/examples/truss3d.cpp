//
// CALFEM C++ example: 3D truss
//
// Three bars connect three fixed supports on the ground to a top node, which is
// loaded with a vertical point load. The example shows the basic CALFEM workflow:
//
//   1. Define nodes, dofs and element topology
//   2. Extract element dofs and coordinates (extract_edof)
//   3. Compute and assemble element stiffness matrices (bar3e, spassem)
//   4. Solve with boundary conditions (spsolveq)
//   5. Compute element normal forces (extract_eldisp, bar3s)
//
// It uses the calfem-python style API, where results are returned:
//
//   auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);
//   auto [a, r] = spsolveq(K, f, bcDofs);
//
// SPDX-License-Identifier: MIT
//

#include <calfem/calfem.h>

#include <iomanip>
#include <iostream>

using namespace calfem;

int main()
{
    // ---- Model definition (all numbering is 1-based, as in CALFEM) ----

    // Node coordinates [m]

    Matrix coords(4, 3);
    coords << 0.0, 0.0, 0.0, // node 1 (support)
        2.0, 0.0, 0.0,       // node 2 (support)
        0.0, 2.0, 0.0,       // node 3 (support)
        0.5, 0.5, 2.0;       // node 4 (loaded top node)

    // Degrees of freedom per node (ux, uy, uz)

    IntMatrix dofs(4, 3);
    dofs << 1, 2, 3, //
        4, 5, 6,     //
        7, 8, 9,     //
        10, 11, 12;

    // Element topology (node numbers)

    IntMatrix topo(3, 2);
    topo << 1, 4, //
        2, 4,     //
        3, 4;

    // Element properties: E [Pa], A [m^2]

    RowVec ep(2);
    ep << 210e9, 1e-4;

    auto nDofs = dofs.size();
    auto nElements = topo.rows();

    // ---- Element dofs and coordinates ----

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    // ---- Assemble global stiffness matrix ----

    TripletList triplets;

    for (Eigen::Index i = 0; i < nElements; i++)
    {
        Matrix Ke = bar3e(ex.row(i), ey.row(i), ez.row(i), ep);
        spassem(edof.row(i), triplets, Ke);
    }

    SpMatrix K(nDofs, nDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());

    // ---- Loads and boundary conditions ----

    ColVec f = ColVec::Zero(nDofs);

    double P = -10e3; // 10 kN downwards at node 4
    f(dofs(3, 2) - 1) += P;

    IntColVec bcDofs(9); // nodes 1-3 fully fixed (prescribed values default to zero)
    bcDofs << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    // ---- Solve ----

    ColVec a, r;

    try
    {
        std::tie(a, r) = spsolveq(K, f, bcDofs);
    }
    catch (const SolverError &e)
    {
        std::cerr << "Solver failed: " << e.what() << "\n";
        return 1;
    }

    // ---- Element normal forces ----

    Matrix ed = extract_eldisp(edof, a);

    std::cout << std::scientific << std::setprecision(4);

    std::cout << "\nDisplacement of node 4 [m]:\n";
    std::cout << "  ux = " << a(dofs(3, 0) - 1) << "\n";
    std::cout << "  uy = " << a(dofs(3, 1) - 1) << "\n";
    std::cout << "  uz = " << a(dofs(3, 2) - 1) << "\n";

    std::cout << "\nNormal forces [N]:\n";

    for (Eigen::Index i = 0; i < nElements; i++)
    {
        ColVec es = bar3s(ex.row(i), ey.row(i), ez.row(i), ep, ed.row(i));
        std::cout << "  bar " << i + 1 << ": N = " << es(0) << "\n";
    }

    // ---- Equilibrium check: support reactions balance the applied load ----

    std::cout << "\nSum of support reactions [N]:\n";
    for (int d = 0; d < 3; d++)
    {
        double sum = 0.0;
        for (Eigen::Index node = 0; node < 3; node++)
            sum += r(dofs(node, d) - 1);
        std::cout << "  R" << "xyz"[d] << " = " << sum << "\n";
    }

    return 0;
}
