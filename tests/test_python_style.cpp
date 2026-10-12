//
// calfem-python style element API (returned values / tuples). Every overload must give
// exactly the same result as the output-argument version it wraps.
//

#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace calfem;

namespace {

struct BarData {
    RowVec ex, ey, ez, ep, ed;
    double eq = -250.0;

    BarData() : ex(2), ey(2), ez(2), ep(2), ed(6)
    {
        ex << 0.0, 2.0;
        ey << 0.0, 1.0;
        ez << 0.0, 0.5;
        ep << 210e9, 1e-3;
        ed << 0.0, 0.0, 0.0, 1e-4, -2e-4, 5e-5;
    }
};

struct BeamData {
    RowVec ex, ey, ez, eo, ep, ed, eq;

    BeamData() : ex(2), ey(2), ez(2), eo(3), ep(6), ed(12), eq(4)
    {
        ex << 0.0, 3.0;
        ey << 0.0, 0.5;
        ez << 0.0, 1.0;
        eo << 0.0, 0.0, 1.0;
        ep << 210e9, 80e9, 1e-2, 2e-5, 3e-5, 1e-5;
        ed << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-4, -2e-3, 5e-4, 1e-4, 2e-4, -3e-4;
        eq << 100.0, -2000.0, 500.0, 10.0;
    }
};

} // namespace

// Return types follow calfem-python: a single value without eq/nep, a tuple with them.

static_assert(std::is_same_v<decltype(bar3e(RowVec(), RowVec(), RowVec(), RowVec())), Matrix>);
static_assert(std::is_same_v<decltype(bar3e(RowVec(), RowVec(), RowVec(), RowVec(), 0.0)), std::tuple<Matrix, ColVec>>);
static_assert(std::is_same_v<decltype(bar3s(RowVec(), RowVec(), RowVec(), RowVec(), RowVec())), ColVec>);
static_assert(std::is_same_v<decltype(bar3s(RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), 0.0)), ColVec>);
static_assert(std::is_same_v<decltype(bar3s(RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), 0.0, 5)),
                             std::tuple<ColVec, ColVec, ColVec>>);
static_assert(std::is_same_v<decltype(beam3e(RowVec(), RowVec(), RowVec(), RowVec(), RowVec())), Matrix>);
static_assert(std::is_same_v<decltype(beam3e(RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), RowVec())),
                             std::tuple<Matrix, ColVec>>);
static_assert(std::is_same_v<decltype(beam3s(RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), RowVec())), Matrix>);
static_assert(
    std::is_same_v<decltype(beam3s(RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), RowVec())), Matrix>);
static_assert(std::is_same_v<decltype(beam3s(RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), RowVec(), 5)),
                             std::tuple<Matrix, Matrix, ColVec>>);

TEST_CASE("bar3e: python style matches output-argument style", "[bar3e][python-style]")
{
    BarData d;

    Matrix KeRef;
    ColVec feRef;
    bar3e(d.ex, d.ey, d.ez, d.ep, d.eq, KeRef, feRef);

    auto [Ke, fe] = bar3e(d.ex, d.ey, d.ez, d.ep, d.eq);
    CHECK(Ke == KeRef);
    CHECK(fe == feRef);

    Matrix KeOnly = bar3e(d.ex, d.ey, d.ez, d.ep);
    CHECK(KeOnly == KeRef);
}

TEST_CASE("bar3s: python style matches output-argument style", "[bar3s][python-style]")
{
    BarData d;
    int nep = 7;

    ColVec esRef, ediRef, eciRef;
    bar3s(d.ex, d.ey, d.ez, d.ep, d.ed, d.eq, nep, esRef, ediRef, eciRef);

    auto [es, edi, eci] = bar3s(d.ex, d.ey, d.ez, d.ep, d.ed, d.eq, nep);
    CHECK(es == esRef);
    CHECK(edi == ediRef);
    CHECK(eci == eciRef);

    // Without nep: section forces at the two element ends
    ColVec es2Ref, edi2Ref, eci2Ref;
    bar3s(d.ex, d.ey, d.ez, d.ep, d.ed, d.eq, 2, es2Ref, edi2Ref, eci2Ref);
    CHECK(bar3s(d.ex, d.ey, d.ez, d.ep, d.ed, d.eq) == es2Ref);

    // Without eq: no distributed load
    bar3s(d.ex, d.ey, d.ez, d.ep, d.ed, 0.0, 2, es2Ref, edi2Ref, eci2Ref);
    CHECK(bar3s(d.ex, d.ey, d.ez, d.ep, d.ed) == es2Ref);
}

TEST_CASE("beam3e: python style matches output-argument style", "[beam3e][python-style]")
{
    BeamData d;

    Matrix KeRef;
    ColVec feRef;
    beam3e(d.ex, d.ey, d.ez, d.eo, d.ep, d.eq, KeRef, feRef);

    auto [Ke, fe] = beam3e(d.ex, d.ey, d.ez, d.eo, d.ep, d.eq);
    CHECK(Ke == KeRef);
    CHECK(fe == feRef);

    Matrix KeOnly = beam3e(d.ex, d.ey, d.ez, d.eo, d.ep);
    CHECK(KeOnly == KeRef);
}

TEST_CASE("beam3s: python style matches output-argument style", "[beam3s][python-style]")
{
    BeamData d;
    int nep = 11;

    Matrix esRef, ediRef;
    ColVec eciRef;
    beam3s(d.ex, d.ey, d.ez, d.eo, d.ep, d.ed, d.eq, nep, esRef, ediRef, eciRef);

    auto [es, edi, eci] = beam3s(d.ex, d.ey, d.ez, d.eo, d.ep, d.ed, d.eq, nep);
    CHECK(es == esRef);
    CHECK(edi == ediRef);
    CHECK(eci == eciRef);

    // Without nep: section forces at the two element ends
    beam3s(d.ex, d.ey, d.ez, d.eo, d.ep, d.ed, d.eq, 2, esRef, ediRef, eciRef);
    Matrix es2 = beam3s(d.ex, d.ey, d.ez, d.eo, d.ep, d.ed, d.eq);
    CHECK(es2 == esRef);
    CHECK(es2.rows() == 2);
    CHECK(es2.cols() == 6);

    // Without eq: no distributed load
    beam3s(d.ex, d.ey, d.ez, d.eo, d.ep, d.ed, RowVec::Zero(4), 2, esRef, ediRef, eciRef);
    CHECK(beam3s(d.ex, d.ey, d.ez, d.eo, d.ep, d.ed) == esRef);
}

// ---- Dof handling, assembly and solvers ----

namespace {

// Three springs in series (k), nodes 1-2-3-4
Matrix springChain(double k)
{
    Matrix K = Matrix::Zero(4, 4);
    Matrix Ke(2, 2);
    Ke << k, -k, -k, k;
    for (int e = 0; e < 3; e++)
    {
        IntRowVec topo(2);
        topo << e + 1, e + 2;
        assem(topo, K, Ke);
    }
    return K;
}

bool messageContains(const SolverError &e, const std::string &text)
{
    return std::string(e.what()).find(text) != std::string::npos;
}

} // namespace

static_assert(std::is_same_v<decltype(extractElDisp(IntMatrix(), ColVec())), Matrix>);
static_assert(std::is_same_v<decltype(extractEdof(IntMatrix(), Matrix(), IntMatrix())),
                             std::tuple<IntMatrix, Matrix, Matrix, Matrix>>);
static_assert(std::is_same_v<decltype(extractAllDofs(0, std::vector<Eigen::Index>())),
                             std::tuple<std::vector<Eigen::Index>, std::vector<Eigen::Index>>>);
static_assert(std::is_same_v<decltype(solveq(Matrix(), ColVec(), IntColVec())), std::tuple<ColVec, ColVec>>);
static_assert(
    std::is_same_v<decltype(solveq(Matrix(), ColVec(), IntColVec(), ColVec())), std::tuple<ColVec, ColVec>>);
static_assert(std::is_same_v<decltype(spsolveq(SpMatrix(), ColVec(), IntColVec())), std::tuple<ColVec, ColVec>>);
static_assert(
    std::is_same_v<decltype(spsolveq(SpMatrix(), ColVec(), IntColVec(), ColVec())), std::tuple<ColVec, ColVec>>);
static_assert(std::is_same_v<decltype(std::declval<SparseSolver &>().solve(ColVec())), std::tuple<ColVec, ColVec>>);
static_assert(
    std::is_same_v<decltype(std::declval<SparseSolver &>().recompute(ColVec())), std::tuple<ColVec, ColVec>>);

TEST_CASE("extractAllDofs/extractElDisp/extractEdof: python style", "[dofs][python-style]")
{
    std::vector<Eigen::Index> freeRef, allRef;
    extractAllDofs(6, {5, 2}, freeRef, allRef);

    auto [freeDofs, allDofs] = extractAllDofs(6, {5, 2});
    CHECK(freeDofs == freeRef);
    CHECK(allDofs == allRef);

    Matrix coords(3, 3);
    coords << 0, 0, 0, 1, 0, 0, 2, 1, 0;
    IntMatrix dofs(3, 2);
    dofs << 1, 2, 3, 4, 5, 6;
    IntMatrix topo(2, 2);
    topo << 1, 2, 2, 3;

    IntMatrix edofRef;
    Matrix exRef, eyRef, ezRef;
    extractEdof(dofs, coords, topo, edofRef, exRef, eyRef, ezRef);

    auto [edof, ex, ey, ez] = extractEdof(dofs, coords, topo);
    CHECK(edof == edofRef);
    CHECK(ex == exRef);
    CHECK(ey == eyRef);
    CHECK(ez == ezRef);

    ColVec a(6);
    a << 10, 20, 30, 40, 50, 60;

    Matrix edRef;
    extractElDisp(edof, a, edRef);
    CHECK(extractElDisp(edof, a) == edRef);
}

TEST_CASE("assem/spassem: stiffness only", "[assem][python-style]")
{
    Matrix Ke(2, 2);
    Ke << 1.0, -2.0, -3.0, 4.0; // deliberately unsymmetric to check orientation
    ColVec fe(2);
    fe << 5.0, 6.0;
    IntRowVec topo(2);
    topo << 3, 1;

    Matrix KRef = Matrix::Zero(3, 3);
    ColVec fRef = ColVec::Zero(3);
    assem(topo, KRef, Ke, fRef, fe);

    Matrix K = Matrix::Zero(3, 3);
    assem(topo, K, Ke);
    CHECK(K == KRef);

    TripletList tripletsRef, triplets;
    ColVec fSp = ColVec::Zero(3);
    spassem(topo, tripletsRef, Ke, fSp, fe);
    spassem(topo, triplets, Ke);

    SpMatrix KSpRef(3, 3), KSp(3, 3);
    KSpRef.setFromTriplets(tripletsRef.begin(), tripletsRef.end());
    KSp.setFromTriplets(triplets.begin(), triplets.end());
    CHECK(Matrix(KSp) == Matrix(KSpRef));
    CHECK(Matrix(KSp) == KRef);
}

TEST_CASE("solveq/spsolveq: python style", "[solveq][spsolveq][python-style]")
{
    double k = 50.0;
    Matrix K = springChain(k);
    SpMatrix KSp = K.sparseView();
    ColVec f = ColVec::Zero(4);
    f(2) = 10.0;

    IntColVec bcDofs(2);
    bcDofs << 4, 1;
    ColVec bcVals(2);
    bcVals << 0.3, 0.0;

    ColVec aRef, rRef;
    REQUIRE(solveq(K, f, bcDofs, bcVals, aRef, rRef));

    auto [a, r] = solveq(K, f, bcDofs, bcVals);
    CHECK(a == aRef);
    CHECK(r == rRef);

    auto [aSp, rSp] = spsolveq(KSp, f, bcDofs, bcVals);
    CHECK(aSp.isApprox(aRef, 1e-12));
    CHECK(rSp.isApprox(rRef, 1e-12));

    // Without bcVals: prescribed values are zero
    REQUIRE(solveq(K, f, bcDofs, ColVec::Zero(2), aRef, rRef));

    auto [a0, r0] = solveq(K, f, bcDofs);
    CHECK(a0 == aRef);
    CHECK(r0 == rRef);

    auto [a0Sp, r0Sp] = spsolveq(KSp, f, bcDofs);
    CHECK(a0Sp.isApprox(aRef, 1e-12));
}

TEST_CASE("solveq/spsolveq: python style throws SolverError", "[solveq][spsolveq][python-style]")
{
    Matrix K = springChain(50.0);
    SpMatrix KSp = K.sparseView();
    ColVec f = ColVec::Zero(4);

    IntColVec outOfRange(1);
    outOfRange << 7;

    try
    {
        solveq(K, f, outOfRange);
        FAIL("solveq did not throw");
    }
    catch (const SolverError &e)
    {
        CHECK(messageContains(e, "outside 1..4"));
    }

    try
    {
        spsolveq(KSp, f, outOfRange);
        FAIL("spsolveq did not throw");
    }
    catch (const SolverError &e)
    {
        CHECK(messageContains(e, "outside 1..4"));
    }

    // Dof 5 has no stiffness: factorization fails
    Matrix K5 = Matrix::Zero(5, 5);
    K5.topLeftCorner(4, 4) = K;
    IntColVec bcDofs(1);
    bcDofs << 1;

    CHECK_THROWS_AS(solveq(K5, ColVec::Zero(5), bcDofs), SolverError);
    CHECK_THROWS_AS(spsolveq(SpMatrix(K5.sparseView()), ColVec::Zero(5), bcDofs), SolverError);

    // SolverError is a std::runtime_error
    CHECK_THROWS_AS(solveq(K, f, outOfRange), std::runtime_error);
}

TEST_CASE("SparseSolver: python style", "[SparseSolver][python-style]")
{
    double k = 50.0;
    SpMatrix K = springChain(k).sparseView();
    ColVec f = ColVec::Zero(4);
    f(1) = 10.0;

    IntColVec bcDofs(2);
    bcDofs << 1, 4;

    ColVec aRef, rRef;
    REQUIRE(spsolveq(K, f, bcDofs, ColVec::Zero(2), aRef, rRef));

    SparseSolver solver;
    REQUIRE(solver.setup(K, bcDofs));

    auto [a, r] = solver.solve(f);
    CHECK(a.isApprox(aRef, 1e-12));
    CHECK(r.isApprox(rRef, 1e-12));

    auto [a2, r2] = solver.recompute(2.0 * f);
    CHECK(a2.isApprox(2.0 * aRef, 1e-12));

    SparseSolver notSetUp;
    CHECK_THROWS_AS(notSetUp.solve(f), SolverError);
    CHECK_THROWS_AS(notSetUp.recompute(f), SolverError);
}

TEST_CASE("extract_eldisp/extract_edof: calfem-python names", "[dofs][python-style]")
{
    Matrix coords(3, 3);
    coords << 0, 0, 0, 1, 0, 0, 2, 1, 0;
    IntMatrix dofs(3, 2);
    dofs << 1, 2, 3, 4, 5, 6;
    IntMatrix topo(2, 2);
    topo << 1, 2, 2, 3;

    IntMatrix edofRef;
    Matrix exRef, eyRef, ezRef;
    extractEdof(dofs, coords, topo, edofRef, exRef, eyRef, ezRef);

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);
    CHECK(edof == edofRef);
    CHECK(ex == exRef);
    CHECK(ey == eyRef);
    CHECK(ez == ezRef);

    IntMatrix edof2;
    Matrix ex2, ey2, ez2;
    extract_edof(dofs, coords, topo, edof2, ex2, ey2, ez2);
    CHECK(edof2 == edofRef);
    CHECK(ex2 == exRef);

    ColVec a(6);
    a << 10, 20, 30, 40, 50, 60;

    Matrix edRef = extractElDisp(edof, a);
    CHECK(extract_eldisp(edof, a) == edRef);

    Matrix ed2;
    extract_eldisp(edof, a, ed2);
    CHECK(ed2 == edRef);
}
