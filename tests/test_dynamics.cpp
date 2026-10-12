//
// Time integration, eigenvalues and condensation: checks not covered by the calfem-python
// reference data (analytical solutions, error handling, deviations from calfem-python).
//

#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace calfem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const double pi = 3.14159265358979323846;

RowVec row(std::initializer_list<double> values)
{
    RowVec v(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for (auto value : values)
        v(i++) = value;
    return v;
}

Matrix springChain(int n, double k)
{
    Matrix K = Matrix::Zero(n, n);
    for (int i = 0; i < n; i++)
    {
        K(i, i) = (i < n - 1) ? 2.0 * k : k;
        if (i > 0)
            K(i, i - 1) = K(i - 1, i) = -k;
    }
    return K;
}

} // namespace

TEST_CASE("step2: undamped single dof oscillator", "[dynamics]")
{
    // m d2a + k a = 0, a(0) = 1: a(t) = cos(w t). The average acceleration method
    // conserves the amplitude exactly and has a small period error for small dt.
    double k = 4.0 * pi * pi, m = 1.0; // w = 2 pi, period 1 s
    Matrix K(1, 1), M(1, 1);
    K << k;
    M << m;

    ColVec a0(1), da0(1);
    a0 << 1.0;
    da0 << 0.0;

    double dt = 0.001;
    auto [model, dof] = step2(K, Matrix(), M, Matrix(), a0, da0, Matrix(), row({dt, 1.0, 0.25, 0.5}), RowVec(),
                              IntColVec());

    REQUIRE(model.a.cols() == 1001); // 0..1000 steps (int(1/0.001) = 1000)
    CHECK(dof.a.size() == 0);

    // After one period the displacement is back at 1
    CHECK_THAT(model.a(0, 1000), WithinAbs(1.0, 1e-4));
    // Half a period: -1
    CHECK_THAT(model.a(0, 500), WithinAbs(-1.0, 1e-4));

    // Energy is conserved by the trapezoidal rule
    for (Eigen::Index j = 0; j < model.a.cols(); j += 100)
    {
        double energy = 0.5 * k * model.a(0, j) * model.a(0, j) + 0.5 * m * model.da(0, j) * model.da(0, j);
        CHECK_THAT(energy, WithinRel(0.5 * k, 1e-10));
    }
}

TEST_CASE("step1: exponential decay", "[dynamics]")
{
    // c da + k a = 0: a(t) = a0 exp(-k t / c). Crank-Nicolson is second order accurate.
    Matrix K(1, 1), C(1, 1);
    K << 2.0;
    C << 1.0;
    ColVec a0(1);
    a0 << 1.0;

    double dt = 0.001;
    auto [model, dof] = step1(K, C, Matrix(), a0, Matrix(), row({dt, 1.0, 0.5}), RowVec(), IntColVec());

    double tEnd = dt * static_cast<double>(model.a.cols() - 1);
    CHECK_THAT(model.a(0, model.a.cols() - 1), WithinRel(std::exp(-2.0 * tEnd), 1e-6));
}

TEST_CASE("step1/step2: dofhist is recorded without times", "[dynamics]")
{
    // Deviation from calfem-python, which leaves dofhist at zero after t = 0 in this case
    Matrix K = springChain(4, 100.0);
    Matrix C = Matrix::Identity(4, 4);
    Matrix f(4, 1);
    f << 0.0, 1.0, 0.0, -2.0;
    ColVec a0 = ColVec::Zero(4);

    IntColVec dofs(2);
    dofs << 2, 4;

    auto [model1, dof1] = step1(K, C, f, a0, Matrix(), row({0.01, 0.5, 0.5}), RowVec(), dofs);

    REQUIRE(dof1.a.rows() == 2);
    REQUIRE(dof1.a.cols() == model1.a.cols());
    CHECK(dof1.a.row(0) == model1.a.row(1));
    CHECK(dof1.a.row(1) == model1.a.row(3));
    CHECK(dof1.da.row(1) == model1.da.row(3));
    CHECK_FALSE(dof1.a.isZero());

    auto [model2, dof2] = step2(K, Matrix(), C, f, a0, a0, Matrix(), row({0.01, 0.5, 0.25, 0.5}), RowVec(), dofs);

    CHECK(dof2.a.row(0) == model2.a.row(1));
    CHECK(dof2.d2a.row(1) == model2.d2a.row(3));
}

TEST_CASE("step1/step2: prescribed values and stored times", "[dynamics]")
{
    Matrix K = springChain(3, 10.0);
    Matrix C = Matrix::Identity(3, 3);
    ColVec a0 = ColVec::Zero(3);

    Matrix bc(1, 2);
    bc << 3, 0.5; // dof 3 held at 0.5

    RowVec times = row({0.0, 0.25, 0.5});

    auto [model, dof] = step1(K, C, Matrix(), a0, bc, row({0.05, 1.0, 1.0}), times, IntColVec());

    REQUIRE(model.a.cols() == 3);
    CHECK(model.a(2, 0) == 0.0); // t = 0 stores a0
    CHECK(model.a(2, 1) == 0.5);
    CHECK(model.a(2, 2) == 0.5);
}

TEST_CASE("step1/step2: invalid input throws SolverError", "[dynamics]")
{
    Matrix K = springChain(2, 10.0);
    ColVec a0 = ColVec::Zero(2);

    CHECK_THROWS_AS(step1(K, Matrix::Zero(2, 2), Matrix(), a0, Matrix(), row({0.1, 1.0, 0.5}), RowVec(), IntColVec()),
                    SolverError);
    CHECK_THROWS_AS(step1(K, Matrix::Identity(2, 2), Matrix(), a0, Matrix(), row({0.1, 1.0}), RowVec(), IntColVec()),
                    SolverError);
    CHECK_THROWS_AS(step1(K, Matrix::Identity(2, 2), Matrix(), a0, Matrix(), row({0.0, 1.0, 0.5}), RowVec(),
                          IntColVec()),
                    SolverError);

    IntColVec badDofs(1);
    badDofs << 3;
    CHECK_THROWS_AS(step2(K, Matrix(), Matrix::Identity(2, 2), Matrix(), a0, a0, Matrix(),
                          row({0.1, 1.0, 0.25, 0.5}), RowVec(), badDofs),
                    SolverError);
    CHECK_THROWS_AS(step2(K, Matrix(), Matrix::Zero(2, 2), Matrix(), a0, a0, Matrix(), row({0.1, 1.0, 0.25, 0.5}),
                          RowVec(), IntColVec()),
                    SolverError);
}

TEST_CASE("eigen: two-dof spring-mass system", "[eigen]")
{
    // Fixed - k - m - k - m: eigenvalues k/m (3 -+ sqrt(5)) / 2
    double k = 3.0, m = 2.0;
    Matrix K = springChain(3, k); // dof 1 is the fixed end
    Matrix M = Matrix::Identity(3, 3) * m;

    IntColVec b(1);
    b << 1;

    auto [L, X] = eigen(K, M, b);

    REQUIRE(L.size() == 2);
    CHECK_THAT(L(0), WithinRel(k / m * (3.0 - std::sqrt(5.0)) / 2.0, 1e-12));
    CHECK_THAT(L(1), WithinRel(k / m * (3.0 + std::sqrt(5.0)) / 2.0, 1e-12));

    // Fixed dof has zero amplitude, modes are mass-normalized
    CHECK(X.row(0).isZero());
    CHECK((X.transpose() * M * X).isApprox(Matrix::Identity(2, 2), 1e-12));

    // K X = M X L on the free dofs
    Matrix Kf = red(K, b), Mf = red(M, b);
    Matrix Xf = X.bottomRows(2);
    CHECK((Kf * Xf - Mf * Xf * L.asDiagonal()).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("eigen/statcon/red: invalid input throws SolverError", "[eigen][statcon][red]")
{
    Matrix K = springChain(3, 1.0);

    IntColVec bad(1);
    bad << 4;

    CHECK_THROWS_AS(red(K, bad), SolverError);
    CHECK_THROWS_AS(statcon(K, ColVec::Zero(3), bad), SolverError);
    CHECK_THROWS_AS(eigen(K, Matrix::Identity(2, 2)), SolverError);
    CHECK_THROWS_AS(eigen(K, -Matrix::Identity(3, 3)), SolverError);
}

TEST_CASE("statcon: condensed system gives the same solution", "[statcon]")
{
    Matrix K = springChain(4, 10.0);
    ColVec f(4);
    f << 1.0, 2.0, 3.0, 4.0;

    IntColVec cd(2);
    cd << 2, 3;

    auto [K1, f1] = statcon(K, f, cd);
    ColVec a1 = K1.ldlt().solve(f1);
    ColVec a = K.ldlt().solve(f);

    CHECK_THAT(a1(0), WithinRel(a(0), 1e-12));
    CHECK_THAT(a1(1), WithinRel(a(3), 1e-12));
}

TEST_CASE("gfunc: interpolation", "[gfunc]")
{
    Matrix G(3, 2);
    G << 0.0, 0.0, //
        1.0, 2.0,  //
        2.0, 0.0;

    auto [ti, g] = gfunc(G, 0.5);

    REQUIRE(ti.size() == 5);
    CHECK(g(1) == 1.0);
    CHECK(g(2) == 2.0);
    CHECK(g(3) == 1.0);
    CHECK(g(4) == 0.0);

    // dt = 0.3 over [0, 2]: points up to 1.8, never beyond the last time
    auto [ti2, g2] = gfunc(G, 0.3);
    CHECK(ti2.size() == 7);
    CHECK(ti2(6) <= 2.0);
}
