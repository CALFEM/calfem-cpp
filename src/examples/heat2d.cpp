//
// CALFEM C++ example: steady-state heat conduction in a plate
//
// A 2 m x 1 m plate with uniform heat generation Q is held at 0 C on the left edge and
// 100 C on the right edge; the top and bottom edges are insulated. The temperature only
// varies in x, with the exact solution
//
//   T(x) = Q / (2 k) x (L - x) + T1 x / L
//
// The plate is meshed with 4-node isoparametric elements (flw2i4e), assembled sparsely and
// solved with spsolveq. The heat flux is computed with flw2i4s.
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

    double Lx = 2.0, Ly = 1.0;
    double k = 1.5;   // thermal conductivity [W/(m C)]
    double Q = 20.0;  // heat generation [W/m^3]
    double T1 = 100.0; // temperature at x = Lx [C]

    int nx = 20, ny = 4;

    // ---- Structured mesh, one dof (temperature) per node ----

    int nnodes = (nx + 1) * (ny + 1);
    auto node = [nx](int i, int j) { return j * (nx + 1) + i + 1; }; // 1-based

    Matrix coords(nnodes, 2);
    IntMatrix dofs(nnodes, 1);
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
        {
            int n = node(i, j) - 1;
            coords(n, 0) = Lx * i / nx;
            coords(n, 1) = Ly * j / ny;
            dofs(n, 0) = n + 1;
        }

    IntMatrix topo(nx * ny, 4);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            topo.row(j * nx + i) << node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1);

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

    // ---- Assembly ----

    RowVec ep(2);
    ep << 1.0, 2; // thickness, 2 x 2 Gauss points
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

    // ---- Boundary conditions: T = 0 on the left edge, T = T1 on the right edge ----

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

    // ---- Compare with the exact solution ----

    auto exact = [&](double x) { return Q / (2.0 * k) * x * (Lx - x) + T1 * x / Lx; };

    double maxError = 0.0;
    for (int n = 0; n < nnodes; n++)
        maxError = std::max(maxError, std::abs(a(n) - exact(coords(n, 0))));

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nTemperature along y = 0 [C]:\n";
    std::cout << "      x        FE        exact\n";
    for (int i = 0; i <= nx; i += 4)
    {
        double x = Lx * i / nx;
        std::cout << "  " << std::setw(6) << x << "  " << std::setw(9) << a(node(i, 0) - 1) << "  " << std::setw(9)
                  << exact(x) << "\n";
    }

    std::cout << std::scientific << std::setprecision(2);
    std::cout << "\nMaximum nodal error: " << maxError << " C\n";

    // ---- Heat flux in the first element (at the left edge) ----

    Matrix ed = extract_eldisp(edof, a);
    auto [es, et, eci] = flw2i4s(ex.row(0), ey.row(0), ep, D, ed.row(0));

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nHeat flux qx at the Gauss points of element 1 [W/m^2]:\n";
    for (Eigen::Index i = 0; i < es.rows(); i++)
    {
        double x = eci(i, 0);
        double qExact = -k * (Q / (2.0 * k) * (Lx - 2.0 * x) + T1 / Lx);
        std::cout << "  x = " << x << ": " << es(i, 0) << " (exact " << qExact << ")\n";
    }

    // Total heat leaving through the fixed edges equals the heat generated
    std::cout << "\nHeat generated:          " << Q * Lx * Ly << " W/m\n";
    std::cout << "Heat through fixed edges: " << -r(bcDofs.array() - 1).sum() << " W/m\n";

    return 0;
}
