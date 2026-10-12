//
// Plane and solid elements, materials: checks not covered by the calfem-python reference data.
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

// Nodal displacements of the linear field u = [1e-3 x + 2e-4 y, -5e-4 x + 3e-3 y] (+ z terms)
// with strains exx = 1e-3, eyy = 3e-3, gxy = 2e-4 - 5e-4 = -3e-4.
RowVec linearDisplacements2D(const RowVec &ex, const RowVec &ey)
{
    RowVec ed(2 * ex.size());
    for (Eigen::Index i = 0; i < ex.size(); i++)
    {
        ed(2 * i) = 1e-3 * ex(i) + 2e-4 * ey(i);
        ed(2 * i + 1) = -5e-4 * ex(i) + 3e-3 * ey(i);
    }
    return ed;
}

const double E = 210e9, nu = 0.3;

} // namespace

TEST_CASE("hooke: plane stress condensation of the 4 x 4 matrix", "[material]")
{
    // A 4 x 4 D with ptype 1 must give the same element as the 3 x 3 plane stress D
    RowVec ex = row({0.0, 1.0, 0.2}), ey = row({0.0, 0.1, 0.9});
    RowVec ep = row({1, 0.1});

    Matrix K3 = plante(ex, ey, ep, hooke(1, E, nu));
    Matrix K4 = plante(ex, ey, ep, hooke(2, E, nu));

    CHECK((K3 - K4).cwiseAbs().maxCoeff() < 1e-10 * K3.cwiseAbs().maxCoeff());
}

TEST_CASE("hooke/effmises: basic properties", "[material]")
{
    for (int ptype : {1, 2, 3, 4})
    {
        Matrix D = hooke(ptype, E, nu);
        CHECK(D.isApprox(D.transpose()));
    }

    // Uniaxial stress: effective stress equals the stress
    Matrix es(1, 6);
    es << 100.0, 0, 0, 0, 0, 0;
    CHECK_THAT(effmises(es, 4)(0), WithinRel(100.0, 1e-14));

    // Pure shear: sqrt(3) tau
    Matrix tau(1, 3);
    tau << 0, 0, 10.0;
    CHECK_THAT(effmises(tau, 1)(0), WithinRel(std::sqrt(3.0) * 10.0, 1e-14));

    CHECK_THROWS_AS(hooke(5, E, nu), std::invalid_argument);
    CHECK_THROWS_AS(effmises(Matrix::Zero(1, 5), 4), std::invalid_argument);
    CHECK_THROWS_AS(effmises(Matrix::Zero(1, 3), 2), std::invalid_argument);
}

TEST_CASE("Plane elements reproduce a constant strain field", "[plane]")
{
    RowVec expected = row({1e-3, 3e-3, -3e-4});

    for (int ptype : {1, 2})
    {
        Matrix D = hooke(ptype, E, nu);
        INFO("ptype " << ptype);

        RowVec ex3 = row({0.1, 1.3, 0.4}), ey3 = row({0.0, 0.2, 1.1});
        auto [es, et] = plants(ex3, ey3, row({double(ptype), 0.1}), D, linearDisplacements2D(ex3, ey3));
        CHECK(es.isApprox(et * D.transpose()));
        if (ptype == 1)
            CHECK((et - expected).cwiseAbs().maxCoeff() < 1e-15);

        RowVec ex4 = row({0.0, 1.2, 1.1, -0.1}), ey4 = row({0.0, 0.1, 0.9, 1.0});
        RowVec ed4 = linearDisplacements2D(ex4, ey4);

        auto [esq, etq] = planqs(ex4, ey4, row({double(ptype), 0.1}), D, ed4);
        auto [esi, eti] = plani4s(ex4, ey4, row({double(ptype), 0.1, 2}), D, ed4);

        // In-plane strains [xx, yy, xy]
        auto inPlane = [](const Matrix &e) {
            return e.cols() == 3 ? Matrix(e) : Matrix(e(Eigen::all, std::vector<Eigen::Index>{0, 1, 3}));
        };

        CHECK((inPlane(etq) - expected).cwiseAbs().maxCoeff() < 1e-14);
        for (Eigen::Index i = 0; i < eti.rows(); i++)
            CHECK((inPlane(eti.row(i)) - expected).cwiseAbs().maxCoeff() < 1e-14);
    }
}

TEST_CASE("plants: several load cases of one element", "[plane]")
{
    // calfem-python only evaluates the first row in this case
    RowVec ex = row({0.0, 1.0, 0.2}), ey = row({0.0, 0.1, 0.9});
    RowVec ep = row({1, 0.1});
    Matrix D = hooke(1, E, nu);

    Matrix ed(3, 6);
    ed.setRandom();
    ed *= 1e-3;

    auto [es, et] = plants(ex, ey, ep, D, ed);
    REQUIRE(es.rows() == 3);

    for (Eigen::Index i = 0; i < 3; i++)
    {
        auto [esi, eti] = plants(ex, ey, ep, D, ed.row(i));
        CHECK(es.row(i) == esi.row(0));
        CHECK(et.row(i) == eti.row(0));
    }

    CHECK_THROWS_AS(plants(Matrix::Zero(2, 3), Matrix::Zero(2, 3), ep, D, ed), std::invalid_argument);
}

TEST_CASE("Plane elements: rigid body motion gives no forces", "[plane]")
{
    RowVec ex = row({0.0, 1.2, 1.1, -0.1}), ey = row({0.0, 0.1, 0.9, 1.0});
    Matrix D = hooke(1, E, nu);

    // Translation and small rotation
    ColVec rigid(8);
    for (int i = 0; i < 4; i++)
    {
        rigid(2 * i) = 0.01 - 0.002 * ey(i);
        rigid(2 * i + 1) = -0.02 + 0.002 * ex(i);
    }

    Matrix Kq = planqe(ex, ey, row({1, 0.1}), D);
    Matrix Ki = plani4e(ex, ey, row({1, 0.1, 2}), D);

    CHECK((Kq * rigid).cwiseAbs().maxCoeff() < 1e-12 * Kq.cwiseAbs().maxCoeff());
    CHECK((Ki * rigid).cwiseAbs().maxCoeff() < 1e-12 * Ki.cwiseAbs().maxCoeff());
}

TEST_CASE("plantf: internal forces equal Ke ed", "[plane]")
{
    RowVec ex = row({0.0, 1.0, 0.2}), ey = row({0.0, 0.1, 0.9});
    RowVec ep = row({1, 0.1});
    Matrix D = hooke(1, E, nu);
    RowVec ed = row({1e-4, -2e-4, 3e-4, 1e-4, -1e-4, 2e-4});

    auto [es, et] = plants(ex, ey, ep, D, ed);
    ColVec ef = plantf(ex, ey, ep, es.row(0));
    ColVec Ked = plante(ex, ey, ep, D) * ed.transpose();

    CHECK((ef - Ked).cwiseAbs().maxCoeff() < 1e-10 * Ked.cwiseAbs().maxCoeff());
}

TEST_CASE("soli8e/soli8s: box element", "[solid]")
{
    RowVec ex = row({0, 2, 2, 0, 0, 2, 2, 0}), ey = row({0, 0, 1, 1, 0, 0, 1, 1});
    RowVec ez = row({0, 0, 0, 0, 3, 3, 3, 3});
    Matrix D = hooke(4, E, nu);

    // 2 x 2 x 2 integration is exact for a box: 3 x 3 x 3 gives the same stiffness
    auto [K2, f2] = soli8e(ex, ey, ez, row({2}), D, row({0.0, 0.0, -10.0}));
    auto [K3, f3] = soli8e(ex, ey, ez, row({3}), D, row({0.0, 0.0, -10.0}));

    CHECK((K2 - K3).cwiseAbs().maxCoeff() < 1e-10 * K2.cwiseAbs().maxCoeff());
    CHECK_THAT(f3.sum(), WithinRel(-10.0 * 6.0, 1e-12)); // body force times volume

    // Rigid translation: no forces
    ColVec translation = ColVec::Zero(24);
    for (int i = 0; i < 8; i++)
        translation(3 * i + 1) = 0.01;
    CHECK((K2 * translation).cwiseAbs().maxCoeff() < 1e-12 * K2.cwiseAbs().maxCoeff());

    // Uniform strain field: exx = 1e-3, all other strains zero
    RowVec ed = RowVec::Zero(24);
    for (int i = 0; i < 8; i++)
        ed(3 * i) = 1e-3 * ex(i);

    auto [et, es, eci] = soli8s(ex, ey, ez, row({2}), D, ed);
    REQUIRE(et.rows() == 8);
    for (Eigen::Index i = 0; i < 8; i++)
    {
        CHECK_THAT(et(i, 0), WithinRel(1e-3, 1e-12));
        CHECK_THAT(es(i, 0), WithinRel(D(0, 0) * 1e-3, 1e-12));
    }
}

TEST_CASE("Continuum elements reject invalid arguments", "[plane][solid]")
{
    RowVec ex = row({0, 1, 1, 0}), ey = row({0, 0, 1, 1});
    Matrix D = hooke(1, E, nu);

    CHECK_THROWS_AS(plante(row({0, 1, 0}), row({0, 0, 1}), row({3, 0.1}), D), std::invalid_argument);
    CHECK_THROWS_AS(plani4e(ex, ey, row({1, 0.1, 4}), D), std::invalid_argument);
    CHECK_THROWS_AS(planqs(ex, ey, row({1, 0.1}), D, RowVec::Zero(6)), std::invalid_argument);
    CHECK_THROWS_AS(plantf(row({0, 1, 0}), row({0, 0, 1}), row({1, 0.1}), RowVec::Zero(5)), std::invalid_argument);
}

TEST_CASE("stress2nodal: two elements sharing an edge", "[dofs]")
{
    // Two 4-node elements, nodes 2 and 5 shared (2 dofs per node)
    IntMatrix edof(2, 8);
    edof << 1, 2, 3, 4, 9, 10, 7, 8, //
        3, 4, 5, 6, 11, 12, 9, 10;

    ColVec eseff(2);
    eseff << 4.0, 8.0;

    Matrix ev = stress2nodal(eseff, edof);

    REQUIRE(ev.rows() == 2);
    REQUIRE(ev.cols() == 4);
    CHECK(ev(0, 0) == 1.0);       // only element 1: 4 / 4
    CHECK(ev(0, 1) == 1.0 + 2.0); // shared: 4 / 4 + 8 / 4
    CHECK(ev(1, 1) == 2.0);       // only element 2
}
