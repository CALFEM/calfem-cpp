#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace calfem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

RowVec row(std::initializer_list<double> values)
{
    RowVec v(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for (auto value : values)
        v(i++) = value;
    return v;
}

} // namespace

TEST_CASE("bar3e: axial bar along x", "[bar3e]")
{
    double E = 210e9, A = 0.01, L = 2.0, q = 1000.0;

    Matrix Ke;
    ColVec fe;
    bar3e(row({0.0, L}), row({0.0, 0.0}), row({0.0, 0.0}), row({E, A}), q, Ke, fe);

    REQUIRE(Ke.rows() == 6);
    REQUIRE(Ke.cols() == 6);

    double k = E * A / L;
    CHECK_THAT(Ke(0, 0), WithinRel(k));
    CHECK_THAT(Ke(0, 3), WithinRel(-k));
    CHECK_THAT(Ke(3, 3), WithinRel(k));
    CHECK_THAT(Ke(1, 1), WithinAbs(0.0, 1e-9));

    CHECK_THAT(fe(0), WithinRel(q * L / 2.0));
    CHECK_THAT(fe(3), WithinRel(q * L / 2.0));
}

TEST_CASE("bar3e: inclined bar is symmetric and rigid-body free", "[bar3e]")
{
    Matrix Ke;
    ColVec fe;
    bar3e(row({1.0, 2.0}), row({0.5, 2.5}), row({-1.0, 1.0}), row({1.0, 1.0}), 0.0, Ke, fe);

    CHECK(Ke.isApprox(Ke.transpose()));

    ColVec translation(6);
    translation << 0.3, -0.7, 1.1, 0.3, -0.7, 1.1;
    CHECK((Ke * translation).norm() < 1e-12);
}

TEST_CASE("bar3s: normal force from axial elongation", "[bar3s]")
{
    double E = 2.0, A = 3.0, L = 4.0, delta = 0.1;
    int n = 5;

    ColVec es(n), edi(n), eci(n);
    bar3s(row({0.0, L}), row({0.0, 0.0}), row({0.0, 0.0}), row({E, A}), row({0, 0, 0, delta, 0, 0}), 0.0, n, es, edi,
          eci);

    for (int i = 0; i < n; i++)
        CHECK_THAT(es(i), WithinRel(E * A * delta / L));

    CHECK_THAT(edi(0), WithinAbs(0.0, 1e-12));
    CHECK_THAT(edi(n - 1), WithinRel(delta));
    CHECK_THAT(eci(n - 1), WithinRel(L));
}

TEST_CASE("beam3e: local stiffness terms for beam along x", "[beam3e]")
{
    double E = 210e9, G = 80e9, A = 0.01, Iy = 2e-5, Iz = 3e-5, Kv = 1e-5, L = 3.0;
    double qy = -500.0;

    Matrix Ke;
    ColVec fe;
    beam3e(row({0.0, L}), row({0.0, 0.0}), row({0.0, 0.0}), row({0.0, 0.0, 1.0}), row({E, G, A, Iy, Iz, Kv}),
           row({0.0, qy, 0.0, 0.0}), Ke, fe);

    REQUIRE(Ke.rows() == 12);
    CHECK(Ke.isApprox(Ke.transpose()));

    CHECK_THAT(Ke(0, 0), WithinRel(E * A / L));
    CHECK_THAT(Ke(1, 1), WithinRel(12.0 * E * Iz / std::pow(L, 3)));
    CHECK_THAT(Ke(2, 2), WithinRel(12.0 * E * Iy / std::pow(L, 3)));
    CHECK_THAT(Ke(3, 3), WithinRel(G * Kv / L));
    CHECK_THAT(Ke(4, 4), WithinRel(4.0 * E * Iy / L));
    CHECK_THAT(Ke(5, 5), WithinRel(4.0 * E * Iz / L));

    CHECK_THAT(fe(1), WithinRel(qy * L / 2.0));
    CHECK_THAT(fe(5), WithinRel(qy * L * L / 12.0));
    CHECK_THAT(fe(11), WithinRel(-qy * L * L / 12.0));
}
