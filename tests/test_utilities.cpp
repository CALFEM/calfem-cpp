//
// Plate element, dof utilities, boundary conditions and loads: checks not covered by the
// calfem-python reference data.
//

#include <calfem/calfem.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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

TEST_CASE("platre: symmetry, rigid body motion and load resultant", "[plate]")
{
    RowVec ex = row({0.0, 2.0, 2.0, 0.0}), ey = row({0.0, 0.0, 1.5, 1.5});
    Matrix D = hooke(1, 210e9, 0.3);

    auto [Ke, fe] = platre(ex, ey, row({0.02}), D, -1000.0);

    CHECK(Ke.isApprox(Ke.transpose()));

    // Rigid translation in w: dofs [w, theta_x, theta_y] per node
    ColVec w = ColVec::Zero(12);
    for (int i = 0; i < 4; i++)
        w(3 * i) = 1.0;
    CHECK((Ke * w).cwiseAbs().maxCoeff() < 1e-10 * Ke.cwiseAbs().maxCoeff());

    // Vertical load resultant: q times area
    double Fz = fe(0) + fe(3) + fe(6) + fe(9);
    CHECK_THAT(Fz, WithinRel(-1000.0 * 3.0, 1e-14));
}

TEST_CASE("apply_force_node: load on the selected dof", "[boundary]")
{
    // calfem-python loads the following dof instead (it indexes f with the 1-based dof)
    IntMatrix dofs = create_dofs(3, 2); // [[1, 2], [3, 4], [5, 6]]
    ColVec f = ColVec::Zero(6);

    apply_force_node(1, dofs, f, 10.0, 2); // node 1 (zero-based), y direction: dof 4
    CHECK(f(3) == 10.0);
    CHECK(f.sum() == 10.0);

    apply_force_node(2, dofs, f, 1.0); // all dofs of node 2: dofs 5 and 6
    CHECK(f(4) == 1.0);
    CHECK(f(5) == 1.0);
}

TEST_CASE("apply_bc: existing values are kept", "[boundary]")
{
    BoundaryDofs bdofs{{10, {1, 2, 3, 4}}};

    IntColVec bc(1);
    bc << 3;
    ColVec bcVal(1);
    bcVal << 0.5;

    auto [bc2, val2] = apply_bc(bdofs, bc, bcVal, 10, 2.0);

    REQUIRE(bc2.size() == 4);
    CHECK(bc2 == (IntColVec(4) << 1, 2, 3, 4).finished());
    CHECK(val2(2) == 0.5); // dof 3 keeps its value
    CHECK(val2(0) == 2.0);

    // y dofs only
    auto [bc3, val3] = apply_bc(bdofs, IntColVec(), ColVec(), 10, 0.0, 2);
    CHECK(bc3 == (IntColVec(2) << 2, 4).finished());
}

TEST_CASE("Boundary functions reject unknown markers and dimensions", "[boundary]")
{
    BoundaryDofs bdofs{{10, {1, 2, 3, 4}}};
    ColVec f = ColVec::Zero(4);
    IntMatrix dofs = create_dofs(2, 2);

    CHECK_THROWS_AS(apply_bc(bdofs, IntColVec(), ColVec(), 99), std::invalid_argument);
    CHECK_THROWS_AS(apply_bc(bdofs, IntColVec(), ColVec(), 10, 0.0, 3), std::invalid_argument);
    CHECK_THROWS_AS(apply_force(bdofs, f, 99, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(apply_force_3d(bdofs, f, 10, 1.0, 4), std::invalid_argument);
    CHECK_THROWS_AS(apply_bc_node(5, dofs, IntColVec(), ColVec()), std::invalid_argument);
    CHECK_THROWS_AS(apply_force_node(0, dofs, f, 1.0, 3), std::invalid_argument);
}

TEST_CASE("coord_extract inverts extract_edof", "[dofs]")
{
    Matrix coords(4, 2);
    coords << 0, 0, 1, 0, 1, 1, 0, 1;
    IntMatrix dofs = create_dofs(4, 2);

    IntMatrix topo(2, 3);
    topo << 1, 2, 3, //
        1, 3, 4;

    auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);
    auto [ex2, ey2, ez2] = coord_extract(edof, coords, dofs);

    CHECK(ex2 == ex);
    CHECK(ey2 == ey);
    CHECK(ez2.isZero());

    IntMatrix badEdof = edof;
    badEdof(0, 0) = 99;
    CHECK_THROWS_AS(coord_extract(badEdof, coords, dofs), std::invalid_argument);
}

TEST_CASE("beam2crd: undeformed beam gives the element axis", "[beam]")
{
    Matrix ex(1, 2), ey(1, 2);
    ex << 1.0, 4.0;
    ey << 2.0, 6.0;

    auto [excd, eycd] = beam2crd(ex, ey, Matrix::Zero(1, 6), 10.0);

    REQUIRE(excd.cols() == 21);
    CHECK_THAT(excd(0, 0), WithinAbs(1.0, 1e-14));
    CHECK_THAT(eycd(0, 20), WithinAbs(6.0, 1e-14));
    CHECK_THAT(excd(0, 10), WithinAbs(2.5, 1e-14));
}
