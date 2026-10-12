//
// CALFEM C++ example: 2D portal frame
//
// Two columns (4 m) and a beam (6 m), clamped at both supports. The frame carries a
// horizontal point load at the top of the left column and a uniformly distributed load on
// the beam. The example uses 2D beam elements (beam2e, beam2s) with dense assembly
// (assem) and solveq, written with the calfem-python style API.
//
// SPDX-License-Identifier: MIT
//

#include <calfem/calfem.h>

#include <iomanip>
#include <iostream>

using namespace calfem;

int main()
{
    // ---- Model ----

    Matrix coords(4, 2);
    coords << 0.0, 0.0, // node 1 (support)
        0.0, 4.0,       // node 2
        6.0, 4.0,       // node 3
        6.0, 0.0;       // node 4 (support)

    IntMatrix dofs(4, 3); // ux, uy, theta per node
    dofs << 1, 2, 3,      //
        4, 5, 6,          //
        7, 8, 9,          //
        10, 11, 12;

    IntMatrix topo(3, 2);
    topo << 1, 2, // left column
        2, 3,     // beam
        3, 4;     // right column

    RowVec ep(3); // E [Pa], A [m^2], I [m^4]
    ep << 210e9, 45.3e-4, 2510e-8;

    // Distributed load [qX, qY] in local element coordinates [N/m]
    Matrix eq = Matrix::Zero(3, 2);
    eq(1, 1) = -10e3; // beam: 10 kN/m downwards

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    // ---- Assembly ----

    Matrix K = Matrix::Zero(12, 12);
    ColVec f = ColVec::Zero(12);

    for (Eigen::Index i = 0; i < topo.rows(); i++)
    {
        auto [Ke, fe] = beam2e(ex.row(i), ey.row(i), ep, eq.row(i));
        assem(edof.row(i), K, Ke, f, fe);
    }

    f(dofs(1, 0) - 1) += 2e3; // 2 kN horizontal load at node 2

    // ---- Solve (both supports clamped) ----

    IntColVec bcDofs(6);
    bcDofs << 1, 2, 3, 10, 11, 12;

    auto [a, r] = solveq(K, f, bcDofs);

    // ---- Results ----

    Matrix ed = extract_eldisp(edof, a);

    std::cout << std::scientific << std::setprecision(4);

    std::cout << "\nNode displacements [m, rad]:\n";
    for (Eigen::Index n = 0; n < dofs.rows(); n++)
        std::cout << "  node " << n + 1 << ": ux = " << a(dofs(n, 0) - 1) << ", uy = " << a(dofs(n, 1) - 1)
                  << ", theta = " << a(dofs(n, 2) - 1) << "\n";

    const char *names[] = {"left column", "beam", "right column"};

    std::cout << "\nSection forces at the element ends [N, N, Nm] (N, V, M):\n";
    for (Eigen::Index i = 0; i < topo.rows(); i++)
    {
        Matrix es = beam2s(ex.row(i), ey.row(i), ep, ed.row(i), eq.row(i));
        std::cout << "  " << names[i] << ":\n";
        std::cout << "    start: " << es(0, 0) << ", " << es(0, 1) << ", " << es(0, 2) << "\n";
        std::cout << "    end:   " << es(1, 0) << ", " << es(1, 1) << ", " << es(1, 2) << "\n";
    }

    // Maximum field moment in the beam, from 21 evaluation points
    auto [es, edi, eci] = beam2s(ex.row(1), ey.row(1), ep, ed.row(1), eq.row(1), 21);
    Eigen::Index iMax;
    es.col(2).maxCoeff(&iMax);
    std::cout << "\nMaximum beam moment: " << es(iMax, 2) << " Nm at x = " << eci(iMax) << " m\n";

    // ---- Equilibrium check ----

    double Rx = r(0) + r(9), Ry = r(1) + r(10);
    std::cout << "\nSum of support reactions: Rx = " << Rx << " N, Ry = " << Ry << " N\n";
    std::cout << "Total applied load:       Fx = " << 2e3 << " N, Fy = " << -6.0 * 10e3 << " N\n";

    return 0;
}
