#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>
#include <vector>

using namespace calfem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Three springs in series (k), nodes 1-2-3-4, as a sparse matrix.
SpMatrix springChain(double k)
{
    TripletList triplets;
    for (int e = 0; e < 3; e++)
    {
        triplets.emplace_back(e, e, k);
        triplets.emplace_back(e, e + 1, -k);
        triplets.emplace_back(e + 1, e, -k);
        triplets.emplace_back(e + 1, e + 1, k);
    }
    SpMatrix K(4, 4);
    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}

// u1 = 0 and u4 = 0.3 prescribed, listed in reverse order. Exact solution: u = 0, 0.1, 0.2, 0.3.
void checkSpringSolution(const ColVec &a, const ColVec &Q, double k)
{
    REQUIRE(a.size() == 4);
    REQUIRE(Q.size() == 4);

    CHECK_THAT(a(0), WithinAbs(0.0, 1e-14));
    CHECK_THAT(a(1), WithinRel(0.1, 1e-12));
    CHECK_THAT(a(2), WithinRel(0.2, 1e-12));
    CHECK_THAT(a(3), WithinRel(0.3, 1e-12));

    CHECK_THAT(Q(0), WithinRel(-0.1 * k, 1e-12));
    CHECK_THAT(Q(1), WithinAbs(0.0, 1e-10));
    CHECK_THAT(Q(2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(Q(3), WithinRel(0.1 * k, 1e-12));
}

struct LogCapture {
    std::vector<std::string> messages;

    LogCapture()
    {
        setLogFunction([this](const std::string &message) { messages.push_back(message); });
    }

    ~LogCapture()
    {
        setLogFunction(nullptr);
    }
};

} // namespace

TEST_CASE("spsolveq: unsorted bcDofs and unsized outputs", "[spsolveq]")
{
    double k = 50.0;
    SpMatrix K = springChain(k);
    ColVec f = ColVec::Zero(4);

    IntColVec bcDofs(2);
    bcDofs << 4, 1;
    ColVec bcVals(2);
    bcVals << 0.3, 0.0;

    ColVec a, Q;
    REQUIRE(spsolveq(K, f, bcDofs, bcVals, a, Q));
    checkSpringSolution(a, Q, k);
}

TEST_CASE("SparseSolver: unsorted bcDofs", "[SparseSolver]")
{
    double k = 50.0;
    SpMatrix K = springChain(k);
    ColVec f = ColVec::Zero(4);

    IntColVec bcDofs(2);
    bcDofs << 4, 1;
    ColVec bcVals(2);
    bcVals << 0.3, 0.0;

    SparseSolver solver;
    REQUIRE(solver.setup(K, bcDofs, bcVals));
    CHECK(solver.nFreeDofs() == 2);

    ColVec a, Q;
    REQUIRE(solver.solve(f, a, Q));
    checkSpringSolution(a, Q, k);
}

TEST_CASE("SparseSolver: does not depend on caller data after setup", "[SparseSolver]")
{
    double k = 50.0;
    SparseSolver solver;

    {
        SpMatrix K = springChain(k);
        IntColVec bcDofs(2);
        bcDofs << 1, 4;
        ColVec bcVals(2);
        bcVals << 0.0, 0.3;

        REQUIRE(solver.setup(K, bcDofs, bcVals));
    } // K, bcDofs and bcVals destroyed here

    ColVec f = ColVec::Zero(4);
    ColVec a, Q;
    REQUIRE(solver.solve(f, a, Q));
    checkSpringSolution(a, Q, k);
}

TEST_CASE("SparseSolver: repeated solves give identical results", "[SparseSolver]")
{
    double k = 50.0;
    SpMatrix K = springChain(k);
    IntColVec bcDofs(2);
    bcDofs << 1, 4;
    ColVec bcVals(2);
    bcVals << 0.0, 0.3;
    ColVec f = ColVec::Zero(4);

    SparseSolver solver;
    REQUIRE(solver.setup(K, bcDofs, bcVals));

    ColVec a1, Q1;
    REQUIRE(solver.solve(f, a1, Q1));

    for (int i = 0; i < 5; i++)
    {
        ColVec a, Q;
        REQUIRE(solver.solve(f, a, Q));
        CHECK(a == a1);
        REQUIRE(solver.recompute(f, a, Q));
        CHECK(a == a1);
    }
}

TEST_CASE("SparseSolver: recompute factorizes when solve was not called", "[SparseSolver]")
{
    double k = 50.0;
    SpMatrix K = springChain(k);
    IntColVec bcDofs(2);
    bcDofs << 1, 4;
    ColVec bcVals(2);
    bcVals << 0.0, 0.3;

    SparseSolver solver;
    REQUIRE(solver.setup(K, bcDofs, bcVals));
    CHECK_FALSE(solver.isFactorized());

    ColVec f = ColVec::Zero(4);
    ColVec a, Q;
    REQUIRE(solver.recompute(f, a, Q));
    CHECK(solver.isFactorized());
    checkSpringSolution(a, Q, k);
}

TEST_CASE("SparseSolver: solve without setup fails", "[SparseSolver]")
{
    SparseSolver solver;
    ColVec f = ColVec::Zero(4);
    ColVec a, Q;

    CHECK_FALSE(solver.solve(f, a, Q));
    CHECK_FALSE(solver.recompute(f, a, Q));
}

TEST_CASE("Solvers reject invalid input", "[solveq][spsolveq][SparseSolver]")
{
    LogCapture capture;

    SpMatrix K = springChain(50.0);
    Matrix Kd = Matrix(K);
    ColVec f = ColVec::Zero(4);
    ColVec a, Q;
    SparseSolver solver;

    SECTION("dof out of range")
    {
        IntColVec bcDofs(2);
        bcDofs << 1, 5;
        ColVec bcVals = ColVec::Zero(2);

        CHECK_FALSE(solveq(Kd, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(spsolveq(K, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(solver.setup(K, bcDofs, bcVals));
    }

    SECTION("zero dof number")
    {
        IntColVec bcDofs(1);
        bcDofs << 0;
        ColVec bcVals = ColVec::Zero(1);

        CHECK_FALSE(solveq(Kd, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(spsolveq(K, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(solver.setup(K, bcDofs, bcVals));
    }

    SECTION("duplicate dof")
    {
        IntColVec bcDofs(2);
        bcDofs << 1, 1;
        ColVec bcVals = ColVec::Zero(2);

        CHECK_FALSE(solveq(Kd, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(spsolveq(K, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(solver.setup(K, bcDofs, bcVals));
    }

    SECTION("bcDofs/bcVals size mismatch")
    {
        IntColVec bcDofs(2);
        bcDofs << 1, 4;
        ColVec bcVals = ColVec::Zero(1);

        CHECK_FALSE(solveq(Kd, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(spsolveq(K, f, bcDofs, bcVals, a, Q));
        CHECK_FALSE(solver.setup(K, bcDofs, bcVals));
    }

    SECTION("load vector size mismatch")
    {
        IntColVec bcDofs(1);
        bcDofs << 1;
        ColVec bcVals = ColVec::Zero(1);
        ColVec fWrong = ColVec::Zero(3);

        CHECK_FALSE(solveq(Kd, fWrong, bcDofs, bcVals, a, Q));
        CHECK_FALSE(spsolveq(K, fWrong, bcDofs, bcVals, a, Q));
        REQUIRE(solver.setup(K, bcDofs, bcVals));
        CHECK_FALSE(solver.solve(fWrong, a, Q));
    }

    SECTION("dof without stiffness (exactly singular K)")
    {
        // Dof 5 is not connected to anything
        SpMatrix K5(5, 5);
        TripletList triplets;
        for (Eigen::Index k = 0; k < K.outerSize(); ++k)
            for (SpMatrix::InnerIterator it(K, k); it; ++it)
                triplets.emplace_back(it.row(), it.col(), it.value());
        K5.setFromTriplets(triplets.begin(), triplets.end());

        ColVec f5 = ColVec::Zero(5);
        IntColVec bcDofs(1);
        bcDofs << 1;
        ColVec bcVals = ColVec::Zero(1);

        CHECK_FALSE(solveq(Matrix(K5), f5, bcDofs, bcVals, a, Q));
        CHECK_FALSE(spsolveq(K5, f5, bcDofs, bcVals, a, Q));
        REQUIRE(solver.setup(K5, bcDofs, bcVals));
        CHECK_FALSE(solver.solve(f5, a, Q));
    }

    CHECK_FALSE(capture.messages.empty());
}

TEST_CASE("Solvers warn about unsupported structures", "[solveq][spsolveq][SparseSolver][log]")
{
    // No supports: K is singular, but rounding leaves a tiny positive pivot so the
    // factorization itself succeeds. The solvers keep the original behaviour (no failure)
    // and report a warning through the log function.

    SpMatrix K = springChain(50.0);
    ColVec f = ColVec::Zero(4);
    IntColVec bcDofs(0);
    ColVec bcVals(0);
    ColVec a, Q;

    auto hasWarning = [](const std::vector<std::string> &messages) {
        for (auto &message : messages)
            if (message.find("nearly singular") != std::string::npos)
                return true;
        return false;
    };

    {
        LogCapture capture;
        solveq(Matrix(K), f, bcDofs, bcVals, a, Q);
        CHECK(hasWarning(capture.messages));
    }
    {
        LogCapture capture;
        spsolveq(K, f, bcDofs, bcVals, a, Q);
        CHECK(hasWarning(capture.messages));
    }
    {
        LogCapture capture;
        SparseSolver solver;
        REQUIRE(solver.setup(K, bcDofs, bcVals));
        solver.solve(f, a, Q);
        CHECK(hasWarning(capture.messages));
    }
}

TEST_CASE("Supported structures give no singularity warning", "[log]")
{
    LogCapture capture;

    SpMatrix K = springChain(50.0);
    ColVec f = ColVec::Zero(4);
    IntColVec bcDofs(2);
    bcDofs << 1, 4;
    ColVec bcVals(2);
    bcVals << 0.0, 0.3;
    ColVec a, Q;

    REQUIRE(spsolveq(K, f, bcDofs, bcVals, a, Q));

    for (auto &message : capture.messages)
        CHECK(message.find("nearly singular") == std::string::npos);
}

TEST_CASE("Solvers handle all dofs prescribed", "[solveq][spsolveq]")
{
    double k = 50.0;
    SpMatrix K = springChain(k);
    ColVec f = ColVec::Zero(4);

    IntColVec bcDofs(4);
    bcDofs << 1, 2, 3, 4;
    ColVec bcVals(4);
    bcVals << 0.0, 0.1, 0.2, 0.3;

    ColVec a, Q;
    REQUIRE(spsolveq(K, f, bcDofs, bcVals, a, Q));
    checkSpringSolution(a, Q, k);

    REQUIRE(solveq(Matrix(K), f, bcDofs, bcVals, a, Q));
    checkSpringSolution(a, Q, k);

    SparseSolver solver;
    REQUIRE(solver.setup(K, bcDofs, bcVals));
    REQUIRE(solver.solve(f, a, Q));
    checkSpringSolution(a, Q, k);
}

TEST_CASE("Log function receives solver diagnostics and can be removed", "[log]")
{
    std::vector<std::string> messages;
    setLogFunction([&messages](const std::string &message) { messages.push_back(message); });

    double k = 50.0;
    SpMatrix K = springChain(k);
    ColVec f = ColVec::Zero(4);
    IntColVec bcDofs(2);
    bcDofs << 1, 4;
    ColVec bcVals(2);
    bcVals << 0.0, 0.3;
    ColVec a, Q;

    REQUIRE(spsolveq(K, f, bcDofs, bcVals, a, Q));
    CHECK(messages.size() == 1);

    setLogFunction(nullptr);

    REQUIRE(spsolveq(K, f, bcDofs, bcVals, a, Q));
    CHECK(messages.size() == 1);
}

TEST_CASE("extractAllDofs: unsorted and duplicate bcDofs", "[dofs]")
{
    std::vector<Eigen::Index> bcDofs{5, 2, 2};
    std::vector<Eigen::Index> freeDofs, allDofs;

    extractAllDofs(6, bcDofs, freeDofs, allDofs);

    CHECK(freeDofs == std::vector<Eigen::Index>{1, 3, 4, 6});
}

TEST_CASE("bar3s: outputs are sized by the routine", "[bar3s]")
{
    RowVec ex(2), ey(2), ez(2), ep(2), ed(6);
    ex << 0.0, 2.0;
    ey << 0.0, 0.0;
    ez << 0.0, 0.0;
    ep << 1.0, 1.0;
    ed << 0.0, 0.0, 0.0, 0.2, 0.0, 0.0;

    ColVec es, edi, eci;
    bar3s(ex, ey, ez, ep, ed, 0.0, 3, es, edi, eci);

    REQUIRE(es.size() == 3);
    REQUIRE(edi.size() == 3);
    REQUIRE(eci.size() == 3);
    CHECK_THAT(es(1), WithinRel(0.1));
    CHECK_THAT(eci(1), WithinRel(1.0));
}

TEST_CASE("extractEdof: 2D coordinates", "[dofs]")
{
    Matrix coords(3, 2);
    coords << 0, 0, 1, 0, 1, 1;

    IntMatrix dofs(3, 2);
    dofs << 1, 2, 3, 4, 5, 6;

    IntMatrix topo(1, 3);
    topo << 1, 2, 3;

    IntMatrix edof;
    Matrix ex, ey, ez;
    extractEdof(dofs, coords, topo, edof, ex, ey, ez);

    CHECK(ex(0, 2) == 1.0);
    CHECK(ey(0, 2) == 1.0);
    CHECK(ez.isZero());
    CHECK(edof.row(0) == (IntRowVec(6) << 1, 2, 3, 4, 5, 6).finished());
}
