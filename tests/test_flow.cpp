//
// Field problem elements: checks not covered by the calfem-python reference data.
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

// Linear field a = 1 + 2 x - 3 y (+ 0.5 z) evaluated at the nodes
RowVec linearField(const RowVec &ex, const RowVec &ey, const RowVec &ez = RowVec())
{
    RowVec a = (1.0 + 2.0 * ex.array() - 3.0 * ey.array()).matrix();
    if (ez.size() > 0)
        a += 0.5 * ez;
    return a;
}

Matrix D2()
{
    Matrix D(2, 2);
    D << 3.0, 0.5, //
        0.5, 2.0;
    return D;
}

void checkConstantGradient(const Matrix &es, const Matrix &et, const Matrix &D, const RowVec &grad)
{
    for (Eigen::Index i = 0; i < et.rows(); i++)
    {
        INFO("point " << i);
        CHECK((et.row(i) - grad).cwiseAbs().maxCoeff() < 1e-12);
        CHECK((es.row(i) + grad * D.transpose()).cwiseAbs().maxCoeff() < 1e-12);
    }
}

} // namespace

TEST_CASE("Field elements reproduce a linear field exactly", "[flow]")
{
    Matrix D = D2();
    RowVec grad2 = row({2.0, -3.0});

    SECTION("flw2ts")
    {
        RowVec ex = row({0.1, 1.3, 0.4}), ey = row({0.0, 0.2, 1.1});
        auto [es, et] = flw2ts(ex, ey, D, linearField(ex, ey));
        checkConstantGradient(es, et, D, grad2);
    }

    SECTION("flw2qs")
    {
        RowVec ex = row({0.0, 1.2, 1.1, -0.1}), ey = row({0.0, 0.1, 0.9, 1.0});
        auto [es, et] = flw2qs(ex, ey, row({1.0}), D, linearField(ex, ey));
        checkConstantGradient(es, et, D, grad2);
    }

    SECTION("flw2i4s")
    {
        RowVec ex = row({0.0, 1.2, 1.1, -0.1}), ey = row({0.0, 0.1, 0.9, 1.0});
        auto [es, et, eci] = flw2i4s(ex, ey, row({1.0, 2}), D, linearField(ex, ey));
        REQUIRE(es.rows() == 4);
        checkConstantGradient(es, et, D, grad2);
    }

    SECTION("flw2i8s")
    {
        RowVec ex = row({0.0, 1.2, 1.1, -0.1, 0.6, 1.15, 0.5, -0.05});
        RowVec ey = row({0.0, 0.1, 0.9, 1.0, 0.05, 0.5, 0.95, 0.5});
        auto [es, et, eci] = flw2i8s(ex, ey, row({1.0, 3}), D, linearField(ex, ey));
        REQUIRE(es.rows() == 9);
        checkConstantGradient(es, et, D, grad2);
    }

    SECTION("flw3i8s")
    {
        RowVec ex = row({0, 1, 1, 0, 0, 1, 1, 0}), ey = row({0, 0, 1, 1, 0, 0, 1, 1});
        RowVec ez = row({0, 0, 0, 0, 1, 1, 1, 1});
        ex(6) = 1.1;
        ez(6) = 0.9;

        Matrix D3 = Matrix::Identity(3, 3) * 2.0;
        for (int ir : {2, 3})
        {
            auto [es, et, eci] = flw3i8s(ex, ey, ez, row({double(ir)}), D3, linearField(ex, ey, ez));
            REQUIRE(es.rows() == ir * ir * ir);
            checkConstantGradient(es, et, D3, row({2.0, -3.0, 0.5}));
        }
    }
}

TEST_CASE("flw3i8e: 27-point rule matches the 8-point rule on a box", "[flow]")
{
    // For a box-shaped element B^T D B is integrated exactly by both rules. calfem-python's
    // ir = 3 rule does not satisfy this (see flow.h).
    RowVec ex = row({0, 2, 2, 0, 0, 2, 2, 0}), ey = row({0, 0, 1, 1, 0, 0, 1, 1});
    RowVec ez = row({0, 0, 0, 0, 3, 3, 3, 3});

    Matrix D(3, 3);
    D << 2.0, 0.1, 0.0, //
        0.1, 1.0, 0.2,  //
        0.0, 0.2, 3.0;

    auto [K2, f2] = flw3i8e(ex, ey, ez, row({2}), D, 5.0);
    auto [K3, f3] = flw3i8e(ex, ey, ez, row({3}), D, 5.0);

    CHECK((K2 - K3).cwiseAbs().maxCoeff() < 1e-12 * K2.cwiseAbs().maxCoeff());
    CHECK_THAT(f3.sum(), WithinRel(5.0 * 6.0, 1e-12)); // source times volume
    CHECK_THAT(f2.sum(), WithinRel(5.0 * 6.0, 1e-12));
}

TEST_CASE("flw2qs: multiple load cases", "[flow]")
{
    // calfem-python fails for more than one row in ed
    RowVec ex = row({0.0, 1.2, 1.1, -0.1}), ey = row({0.0, 0.1, 0.9, 1.0});
    RowVec ep = row({0.5});
    Matrix D = D2();

    Matrix ed(2, 4);
    ed << 1.0, 2.0, -1.0, 0.5, //
        0.3, -0.2, 0.7, 1.1;

    auto [es, et] = flw2qs(ex, ey, ep, D, ed, 3.0);

    for (Eigen::Index j = 0; j < 2; j++)
    {
        auto [esj, etj] = flw2qs(ex, ey, ep, D, ed.row(j), 3.0);
        CHECK(es.row(j) == esj.row(0));
        CHECK(et.row(j) == etj.row(0));
    }
}

TEST_CASE("Field elements: source vector integrates to source times area", "[flow]")
{
    RowVec ex = row({0.0, 2.0, 2.2, -0.1}), ey = row({0.0, 0.1, 1.5, 1.4});
    double t = 0.5, Q = 4.0;

    // Shoelace area
    double area = 0.0;
    for (int i = 0; i < 4; i++)
        area += ex(i) * ey((i + 1) % 4) - ex((i + 1) % 4) * ey(i);
    area *= 0.5;

    auto [Kq, fq] = flw2qe(ex, ey, row({t}), D2(), Q);
    auto [Ki, fi] = flw2i4e(ex, ey, row({t, 2}), D2(), Q);

    CHECK_THAT(fq.sum(), WithinRel(Q * t * area, 1e-12));
    CHECK_THAT(fi.sum(), WithinRel(Q * t * area, 1e-12));

    // Constant field: no flux
    CHECK((Kq * ColVec::Ones(4)).cwiseAbs().maxCoeff() < 1e-12);
    CHECK((Ki * ColVec::Ones(4)).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("Field elements reject unsupported integration rules", "[flow]")
{
    RowVec ex = row({0, 1, 1, 0}), ey = row({0, 0, 1, 1});

    CHECK_THROWS_AS(flw2i4e(ex, ey, row({1.0, 4}), D2()), std::invalid_argument);
    CHECK_THROWS_AS(flw2i8e(row({0, 1, 1, 0, 0.5, 1, 0.5, 0}), row({0, 0, 1, 1, 0, 0.5, 1, 0.5}), row({1.0, 0}), D2()),
                    std::invalid_argument);
    CHECK_THROWS_AS(flw3i8e(row({0, 1, 1, 0, 0, 1, 1, 0}), row({0, 0, 1, 1, 0, 0, 1, 1}), row({0, 0, 0, 0, 1, 1, 1, 1}),
                            row({4}), Matrix::Identity(3, 3)),
                    std::invalid_argument);
}
