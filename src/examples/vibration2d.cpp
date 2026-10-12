//
// CALFEM C++ example: vibration of a cantilever beam
//
// A 3 m steel cantilever is modelled with 10 dynamic beam elements (beam2de). The natural
// frequencies from eigen are compared with the analytical Euler-Bernoulli values, and the
// response to a short load pulse at the free end (described with gfunc) is computed with
// Newmark time integration (step2).
//
// SPDX-License-Identifier: MIT
//

#include <calfem/calfem.h>

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace calfem;

int main()
{
    const double pi = 3.14159265358979323846;

    // ---- Model: nodes along x, 3 dofs per node (ux, uy, theta) ----

    int nel = 10;
    double L = 3.0;
    double E = 210e9, A = 45.3e-4, I = 2510e-8, m = 35.6; // m: mass per unit length [kg/m]

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

    // ---- Assemble stiffness and mass matrices ----

    Matrix K = Matrix::Zero(ndof, ndof);
    Matrix M = Matrix::Zero(ndof, ndof);

    for (int e = 0; e < nel; e++)
    {
        auto [Ke, Me, Ce] = beam2de(ex.row(e), ey.row(e), ep);
        assem(edof.row(e), K, Ke);
        assem(edof.row(e), M, Me);
    }

    IntColVec bcDofs(3); // clamped at x = 0
    bcDofs << 1, 2, 3;

    // ---- Natural frequencies ----

    auto [lambda, X] = eigen(K, M, bcDofs);

    // Euler-Bernoulli cantilever: w_n = (beta_n L)^2 sqrt(EI / (m L^4))
    const double betaL[] = {1.875104, 4.694091, 7.854757};

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nBending frequencies [Hz]     FE      analytical\n";

    int mode = 0;
    for (Eigen::Index i = 0; (i < lambda.size()) && (mode < 3); i++)
    {
        // Skip axial modes (no transverse displacement at the free end)
        if (std::abs(X(ndof - 2, i)) < 1e-12)
            continue;

        double fe = std::sqrt(lambda(i)) / (2.0 * pi);
        double fa = betaL[mode] * betaL[mode] * std::sqrt(E * I / (m * std::pow(L, 4))) / (2.0 * pi);
        std::cout << "  mode " << mode + 1 << "              " << std::setw(9) << fe << "  " << std::setw(9) << fa
                  << "\n";
        mode++;
    }

    // ---- Response to a load pulse at the free end ----

    double dt = 0.0005, T = 0.2;

    Matrix G(4, 2); // triangular 1 kN pulse during the first 10 ms
    G << 0.0, 0.0, //
        0.005, 1e3, //
        0.01, 0.0, //
        T, 0.0;

    auto [t, g] = gfunc(G, dt);

    Matrix f = Matrix::Zero(ndof, t.size());
    f.row(ndof - 2) = g.transpose(); // uy at the free end

    Matrix bc(3, 2);
    bc << 1, 0.0, //
        2, 0.0,   //
        3, 0.0;

    IntColVec tip(1);
    tip << ndof - 1; // uy at the free end (1-based)

    ColVec a0 = ColVec::Zero(ndof);
    RowVec ip(4);
    ip << dt, T, 0.25, 0.5; // average acceleration

    auto [modelhist, dofhist] = step2(K, Matrix(), M, f, a0, a0, bc, ip, RowVec(), tip);

    Eigen::Index iMax;
    double uMax = dofhist.a.row(0).cwiseAbs().maxCoeff(&iMax);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nTip response to the load pulse:\n";
    std::cout << "  maximum |uy| = " << uMax << " m at t = " << std::fixed << iMax * dt << " s\n";
    std::cout << "  steps: " << dofhist.a.cols() - 1 << ", dt = " << dt << " s\n";

    return 0;
}
