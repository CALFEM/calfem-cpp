#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace calfem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Cantilever of length L along x, built from nEl beam3e elements, clamped at x = 0
// and loaded with a point load P in the y-direction at the free end.
struct Cantilever {
    double E = 210e9, G = 80e9, A = 0.01, Iy = 2e-5, Iz = 3e-5, Kv = 1e-5;
    double L = 4.0, P = -1000.0;
    int nEl;

    Eigen::Index nDofs = 0;
    SpMatrix K;
    ColVec f;
    IntColVec bcDofs;
    ColVec bcVals;

    explicit Cantilever(int elements = 4) : nEl(elements)
    {
        auto nNodes = nEl + 1;
        nDofs = 6 * nNodes;

        RowVec ep(6), eo(3), eq(4);
        ep << E, G, A, Iy, Iz, Kv;
        eo << 0.0, 0.0, 1.0;
        eq.setZero();

        TripletList triplets;
        f = ColVec::Zero(nDofs);

        double le = L / nEl;

        for (int e = 0; e < nEl; e++)
        {
            RowVec ex(2), ey(2), ez(2);
            ex << e * le, (e + 1) * le;
            ey.setZero();
            ez.setZero();

            Matrix Ke;
            ColVec fe;
            beam3e(ex, ey, ez, eo, ep, eq, Ke, fe);

            IntRowVec topo(12);
            for (int i = 0; i < 12; i++)
                topo(i) = 6 * e + i + 1;

            spassem(topo, triplets, Ke, f, fe);
        }

        K.resize(nDofs, nDofs);
        K.setFromTriplets(triplets.begin(), triplets.end());

        f(6 * nEl + 1) += P;

        bcDofs.resize(6);
        bcDofs << 1, 2, 3, 4, 5, 6;
        bcVals = ColVec::Zero(6);
    }

    double tipDeflection() const
    {
        return P * std::pow(L, 3) / (3.0 * E * Iz);
    }
};

} // namespace

TEST_CASE("spsolveq: cantilever tip deflection", "[spsolveq]")
{
    Cantilever c;

    ColVec a = ColVec::Zero(c.nDofs);
    ColVec Q;

    REQUIRE(spsolveq(c.K, c.f, c.bcDofs, c.bcVals, a, Q));

    CHECK_THAT(a(6 * c.nEl + 1), WithinRel(c.tipDeflection(), 1e-10));

    // Support reactions balance the applied load
    CHECK_THAT(Q(1), WithinRel(-c.P, 1e-10));
    CHECK_THAT(Q(5), WithinRel(-c.P * c.L, 1e-10));
}

TEST_CASE("solveq: cantilever tip deflection matches spsolveq", "[solveq]")
{
    Cantilever c;

    Matrix K = Matrix(c.K);
    ColVec a, Q;

    REQUIRE(solveq(K, c.f, c.bcDofs, c.bcVals, a, Q));

    REQUIRE(a.size() == c.nDofs);
    CHECK_THAT(a(6 * c.nEl + 1), WithinRel(c.tipDeflection(), 1e-10));
    CHECK_THAT(Q(1), WithinRel(-c.P, 1e-10));
    CHECK_THAT(Q(5), WithinRel(-c.P * c.L, 1e-10));

    ColVec aSp = ColVec::Zero(c.nDofs);
    ColVec QSp;
    REQUIRE(spsolveq(c.K, c.f, c.bcDofs, c.bcVals, aSp, QSp));
    CHECK(a.isApprox(aSp, 1e-12));
}

TEST_CASE("solveq: prescribed displacements given in arbitrary order", "[solveq]")
{
    // Two springs in series (k), nodes 1-2-3. u3 = 0.3 is prescribed and u1 = 0,
    // listed in reverse order to check that bcVals follows bcDofs.
    double k = 100.0;

    Matrix K(3, 3);
    K << k, -k, 0.0, //
        -k, 2.0 * k, -k, //
        0.0, -k, k;

    ColVec f = ColVec::Zero(3);

    IntColVec bcDofs(2);
    bcDofs << 3, 1;
    ColVec bcVals(2);
    bcVals << 0.3, 0.0;

    ColVec a, Q;
    REQUIRE(solveq(K, f, bcDofs, bcVals, a, Q));

    CHECK_THAT(a(0), WithinAbs(0.0, 1e-14));
    CHECK_THAT(a(1), WithinRel(0.15, 1e-12));
    CHECK_THAT(a(2), WithinRel(0.3, 1e-12));

    CHECK_THAT(Q(0), WithinRel(-k * 0.15, 1e-12));
    CHECK_THAT(Q(1), WithinAbs(0.0, 1e-12));
    CHECK_THAT(Q(2), WithinRel(k * 0.15, 1e-12));
}

TEST_CASE("solveq: dense assem workflow matches spassem/spsolveq", "[solveq][assem]")
{
    // Three-bar 3D truss (same model as the truss3d example)
    Matrix ex(3, 2), ey(3, 2), ez(3, 2);
    ex << 0.0, 0.5, 2.0, 0.5, 0.0, 0.5;
    ey << 0.0, 0.5, 0.0, 0.5, 2.0, 0.5;
    ez << 0.0, 2.0, 0.0, 2.0, 0.0, 2.0;

    IntMatrix edof(3, 6);
    edof << 1, 2, 3, 10, 11, 12, //
        4, 5, 6, 10, 11, 12,     //
        7, 8, 9, 10, 11, 12;

    RowVec ep(2);
    ep << 210e9, 1e-4;

    Matrix K = Matrix::Zero(12, 12);
    ColVec f = ColVec::Zero(12);
    TripletList triplets;
    ColVec fSp = ColVec::Zero(12);

    for (Eigen::Index i = 0; i < 3; i++)
    {
        Matrix Ke;
        ColVec fe;
        bar3e(ex.row(i), ey.row(i), ez.row(i), ep, 0.0, Ke, fe);

        IntRowVec topo = edof.row(i);
        assem(topo, K, Ke, f, fe);
        spassem(topo, triplets, Ke, fSp, fe);
    }

    SpMatrix KSp(12, 12);
    KSp.setFromTriplets(triplets.begin(), triplets.end());

    CHECK(K.isApprox(Matrix(KSp)));

    f(11) = -10e3;
    fSp(11) = -10e3;

    IntColVec bcDofs(9);
    bcDofs << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    ColVec bcVals = ColVec::Zero(9);

    ColVec a, Q;
    REQUIRE(solveq(K, f, bcDofs, bcVals, a, Q));

    ColVec aSp = ColVec::Zero(12);
    ColVec QSp;
    REQUIRE(spsolveq(KSp, fSp, bcDofs, bcVals, aSp, QSp));

    CHECK(a.isApprox(aSp, 1e-12));
    CHECK(Q.isApprox(QSp, 1e-10));
}

TEST_CASE("SparseSolver: setup/solve/recompute match spsolveq", "[SparseSolver]")
{
    Cantilever c;

    SparseSolver solver;
    REQUIRE(solver.setup(c.K, c.bcDofs, c.bcVals));

    ColVec a, Q;
    REQUIRE(solver.solve(c.f, a, Q));
    CHECK_THAT(a(6 * c.nEl + 1), WithinRel(c.tipDeflection(), 1e-10));

    ColVec f2 = 2.0 * c.f;
    ColVec a2, Q2;
    REQUIRE(solver.recompute(f2, a2, Q2));
    CHECK_THAT(a2(6 * c.nEl + 1), WithinRel(2.0 * c.tipDeflection(), 1e-10));
}

TEST_CASE("beam3s: section forces along a cantilever element", "[beam3s]")
{
    Cantilever c(1);

    ColVec a = ColVec::Zero(c.nDofs);
    ColVec Q;
    REQUIRE(spsolveq(c.K, c.f, c.bcDofs, c.bcVals, a, Q));

    RowVec ep(6), eo(3), eq(4);
    ep << c.E, c.G, c.A, c.Iy, c.Iz, c.Kv;
    eo << 0.0, 0.0, 1.0;
    eq.setZero();

    RowVec ex(2), ey(2), ez(2);
    ex << 0.0, c.L;
    ey.setZero();
    ez.setZero();

    RowVec ed = a.transpose();

    int n = 5;
    Matrix es, edi;
    ColVec eci;
    beam3s(ex, ey, ez, eo, ep, ed, eq, n, es, edi, eci);

    REQUIRE(es.rows() == n);
    REQUIRE(es.cols() == 6);

    for (int i = 0; i < n; i++)
        CHECK_THAT(std::abs(es(i, 1)), WithinRel(std::abs(c.P), 1e-8));

    CHECK_THAT(std::abs(es(0, 5)), WithinRel(std::abs(c.P * c.L), 1e-8));
    CHECK_THAT(es(n - 1, 5), WithinAbs(0.0, 1e-6));
    CHECK_THAT(edi(n - 1, 1), WithinRel(c.tipDeflection(), 1e-8));
}

TEST_CASE("extractAllDofs and toZeroIndex", "[dofs]")
{
    std::vector<Eigen::Index> bcDofs{1, 4};
    std::vector<Eigen::Index> freeDofs, allDofs;

    extractAllDofs(5, bcDofs, freeDofs, allDofs);

    CHECK(allDofs == std::vector<Eigen::Index>{1, 2, 3, 4, 5});
    CHECK(freeDofs == std::vector<Eigen::Index>{2, 3, 5});
    CHECK(toZeroIndex(freeDofs) == std::vector<Eigen::Index>{1, 2, 4});
}

TEST_CASE("extractEdof and extractElDisp", "[dofs]")
{
    Matrix coords(3, 3);
    coords << 0, 0, 0, 1, 0, 0, 2, 1, 0;

    IntMatrix dofs(3, 2);
    dofs << 1, 2, 3, 4, 5, 6;

    IntMatrix topo(2, 2);
    topo << 1, 2, 2, 3;

    IntMatrix edof;
    Matrix ex, ey, ez;
    extractEdof(dofs, coords, topo, edof, ex, ey, ez);

    IntMatrix expectedEdof(2, 4);
    expectedEdof << 1, 2, 3, 4, 3, 4, 5, 6;
    CHECK(edof == expectedEdof);
    CHECK(ex(1, 1) == 2.0);
    CHECK(ey(1, 1) == 1.0);

    ColVec a(6);
    a << 10, 20, 30, 40, 50, 60;

    Matrix ed;
    extractElDisp(edof, a, ed);

    CHECK(ed(1, 0) == 30.0);
    CHECK(ed(1, 3) == 60.0);
}
