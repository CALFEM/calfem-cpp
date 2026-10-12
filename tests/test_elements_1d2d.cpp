//
// Checks for the 1D/2D element routines that are not covered by the calfem-python
// reference data (argument validation, overload consistency, physical sanity).
//

#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>

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

TEST_CASE("Section routines reject nep < 2", "[elements]")
{
    RowVec ex = row({0.0, 2.0}), ey = row({0.0, 0.0});

    CHECK_THROWS_AS(bar1s(ex, row({1.0, 1.0}), row({0.0, 1.0}), 0.0, 1), std::invalid_argument);
    CHECK_THROWS_AS(beam1s(ex, row({1.0, 1.0}), row({0.0, 0.0, 0.0, 0.0}), 0.0, 0), std::invalid_argument);
    CHECK_THROWS_AS(beam2s(ex, ey, row({1.0, 1.0, 1.0}), RowVec::Zero(6), row({0.0, 0.0}), 1), std::invalid_argument);
}

TEST_CASE("beam2de/beam2ds reject invalid ep", "[elements]")
{
    RowVec ex = row({0.0, 2.0}), ey = row({0.0, 0.0});

    CHECK_THROWS_AS(beam2de(ex, ey, row({1.0, 1.0, 1.0})), std::invalid_argument);
    CHECK_THROWS_AS(beam2ds(ex, ey, row({1.0, 1.0, 1.0, 1.0, 1.0}), Matrix::Zero(1, 6), Matrix::Zero(1, 6),
                            Matrix::Zero(1, 6)),
                    std::invalid_argument);
}

TEST_CASE("Short overloads equal the long ones with zero load", "[elements]")
{
    RowVec ex = row({0.3, 2.1}), ey = row({-0.2, 1.4});
    RowVec ep = row({210e9, 1e-2, 2e-5});
    RowVec ed = row({1e-4, -2e-4, 3e-4, 2e-4, 1e-4, -1e-4});

    CHECK(beam2e(ex, ey, ep) == std::get<0>(beam2e(ex, ey, ep, row({0.0, 0.0}))));
    CHECK(beam2s(ex, ey, ep, ed) == std::get<0>(beam2s(ex, ey, ep, ed, row({0.0, 0.0}), 2)));
    CHECK(bar2s(ex, ey, row({210e9, 1e-2}), row({0.0, 0.0, 1e-4, 0.0})) ==
          std::get<0>(bar2s(ex, ey, row({210e9, 1e-2}), row({0.0, 0.0, 1e-4, 0.0}), 0.0, 2)));
    CHECK(beam2gxe(ex, ey, ep, 0.0) == std::get<0>(beam2gxe(ex, ey, ep, 0.0, 0.0)));
}

TEST_CASE("beam2gxe approaches beam2ge for small axial force", "[elements]")
{
    RowVec ex = row({0.0, 3.0}), ey = row({0.0, 0.0});
    RowVec ep = row({210e9, 1e-2, 2e-5});
    double QX = -1e3; // small compared with the buckling load

    Matrix Kge = beam2ge(ex, ey, ep, QX);
    Matrix Kgxe = beam2gxe(ex, ey, ep, QX);

    CHECK((Kge - Kgxe).cwiseAbs().maxCoeff() < 1e-6 * Kge.cwiseAbs().maxCoeff());
}

TEST_CASE("beam1e: cantilever tip deflection", "[elements]")
{
    double E = 210e9, I = 2e-5, L = 3.0, P = -1000.0;

    Matrix Ke = beam1e(row({0.0, L}), row({E, I}));

    // Clamp node 1 (dofs 1, 2), load node 2
    Matrix Kff = Ke.bottomRightCorner(2, 2);
    ColVec f(2);
    f << P, 0.0;
    ColVec a = Kff.ldlt().solve(f);

    CHECK_THAT(a(0), WithinRel(P * L * L * L / (3.0 * E * I), 1e-12));
}

TEST_CASE("bar2gs returns the axial force", "[elements]")
{
    RowVec ex = row({0.0, 3.0}), ey = row({0.0, 4.0}); // L = 5
    RowVec ep = row({1.0, 2.0});
    RowVec ed = row({0.0, 0.0, 0.3, 0.4}); // elongation 0.5 along the bar

    auto [es, QX] = bar2gs(ex, ey, ep, ed);

    CHECK_THAT(QX, WithinRel(2.0 * 0.5 / 5.0, 1e-12));
    CHECK_THAT(es(1), WithinRel(QX, 1e-12));
}
