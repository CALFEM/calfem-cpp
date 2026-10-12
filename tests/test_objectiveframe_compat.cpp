//
// ObjectiveFrame compatibility tests.
//
// ObjectiveFrame (src/ofsolve/beam_solver.cpp) uses the type aliases, bar3e, bar3s, beam3e,
// beam3s, spassem and SparseSolver (setup/solve/recompute, held by value and reused). These
// tests pin the signatures and check the results against a frozen copy of the original code
// (tests/reference, namespace calfem_ref), using the same call patterns as ObjectiveFrame.
//

#include <calfem/calfem.h>
#include <calfem_ref.h>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <random>
#include <string>
#include <type_traits>

using namespace calfem;

// ---- Signatures and types ----

static_assert(std::is_same_v<RowVec, Eigen::Matrix<double, 1, Eigen::Dynamic>>);
static_assert(std::is_same_v<ColVec, Eigen::Matrix<double, Eigen::Dynamic, 1>>);
static_assert(std::is_same_v<Matrix, Eigen::MatrixXd>);
static_assert(std::is_same_v<SpMatrix, Eigen::SparseMatrix<double>>);
static_assert(std::is_same_v<IntRowVec, Eigen::Matrix<int, 1, Eigen::Dynamic>>);
static_assert(std::is_same_v<IntColVec, Eigen::Matrix<int, Eigen::Dynamic, 1>>);
static_assert(std::is_same_v<TripletList, std::vector<Eigen::Triplet<double, std::ptrdiff_t>>>);
static_assert(std::is_same_v<SparseSolverPtr, std::shared_ptr<SparseSolver>>);
static_assert(std::is_default_constructible_v<SparseSolver>);

TEST_CASE("ObjectiveFrame: function signatures", "[compat]")
{
    // Assigning to explicitly typed function pointers fails to compile if a signature changes,
    // while still allowing new overloads to be added.

    void (*pBar3e)(const RowVec &, const RowVec &, const RowVec &, const RowVec &, double, Matrix &, ColVec &) =
        &calfem::bar3e;
    void (*pBar3s)(const RowVec &, const RowVec &, const RowVec &, const RowVec &, const RowVec &, double, int,
                   ColVec &, ColVec &, ColVec &) = &calfem::bar3s;
    void (*pBeam3e)(const RowVec &, const RowVec &, const RowVec &, const RowVec &, const RowVec &, const RowVec &,
                    Matrix &, ColVec &) = &calfem::beam3e;
    void (*pBeam3s)(const RowVec &, const RowVec &, const RowVec &, const RowVec &, const RowVec &, const RowVec &,
                    const RowVec &, const int, Matrix &, Matrix &, ColVec &) = &calfem::beam3s;
    void (*pSpassem)(const IntRowVec &, TripletList &, const Matrix &, ColVec &, const ColVec &) = &calfem::spassem;
    void (*pExtractAllDofs)(Eigen::Index, const std::vector<Eigen::Index> &, std::vector<Eigen::Index> &,
                            std::vector<Eigen::Index> &) = &calfem::extractAllDofs;
    std::vector<Eigen::Index> (*pToZeroIndex)(const std::vector<Eigen::Index> &) = &calfem::toZeroIndex;
    bool (*pSpsolveq)(const SpMatrix &, const ColVec &, const IntColVec &, const ColVec &, ColVec &, ColVec &) =
        &calfem::spsolveq;

    bool (SparseSolver::*pSetup)(const SpMatrix &, const IntColVec &, const ColVec &) = &SparseSolver::setup;
    bool (SparseSolver::*pSolve)(const ColVec &, ColVec &, ColVec &) = &SparseSolver::solve;
    bool (SparseSolver::*pRecompute)(const ColVec &, ColVec &, ColVec &) = &SparseSolver::recompute;
    SparseSolverPtr (*pCreate)() = &SparseSolver::create;

    CHECK(pBar3e != nullptr);
    CHECK(pBar3s != nullptr);
    CHECK(pBeam3e != nullptr);
    CHECK(pBeam3s != nullptr);
    CHECK(pSpassem != nullptr);
    CHECK(pExtractAllDofs != nullptr);
    CHECK(pToZeroIndex != nullptr);
    CHECK(pSpsolveq != nullptr);
    CHECK(pSetup != nullptr);
    CHECK(pSolve != nullptr);
    CHECK(pRecompute != nullptr);
    CHECK(pCreate != nullptr);
}

// ---- Element routines against the reference ----

namespace {

RowVec randomRow(std::mt19937 &rng, Eigen::Index n, double lo, double hi)
{
    std::uniform_real_distribution<double> dist(lo, hi);
    RowVec v(n);
    for (Eigen::Index i = 0; i < n; i++)
        v(i) = dist(rng);
    return v;
}

} // namespace

TEST_CASE("ObjectiveFrame: bar3e/bar3s match reference", "[compat]")
{
    std::mt19937 rng(1234);

    for (int test = 0; test < 50; test++)
    {
        RowVec Ex = randomRow(rng, 2, -5.0, 5.0);
        RowVec Ey = randomRow(rng, 2, -5.0, 5.0);
        RowVec Ez = randomRow(rng, 2, -5.0, 5.0);
        RowVec Ep(2);
        Ep << 210e9, 1e-3 * (1 + test);
        double eq = randomRow(rng, 1, -1e3, 1e3)(0);

        // Pre-sized outputs, as in ObjectiveFrame
        Matrix Ke_b(6, 6), Ke_ref(6, 6);
        ColVec fe_b(6), fe_ref(6);

        bar3e(Ex, Ey, Ez, Ep, eq, Ke_b, fe_b);
        calfem_ref::bar3e(Ex, Ey, Ez, Ep, eq, Ke_ref, fe_ref);

        CHECK(Ke_b == Ke_ref);
        CHECK(fe_b == fe_ref);

        RowVec Ed_b = randomRow(rng, 6, -1e-3, 1e-3);
        int n = 2 + test % 9;

        ColVec Es_b(n), Edi_b(n), Eci_b(n);
        ColVec Es_ref(n), Edi_ref(n), Eci_ref(n);

        bar3s(Ex, Ey, Ez, Ep, Ed_b, eq, n, Es_b, Edi_b, Eci_b);
        calfem_ref::bar3s(Ex, Ey, Ez, Ep, Ed_b, eq, n, Es_ref, Edi_ref, Eci_ref);

        CHECK(Es_b == Es_ref);
        CHECK(Edi_b == Edi_ref);
        CHECK(Eci_b == Eci_ref);
    }
}

TEST_CASE("ObjectiveFrame: beam3e/beam3s match reference", "[compat]")
{
    std::mt19937 rng(5678);

    for (int test = 0; test < 50; test++)
    {
        RowVec Ex = randomRow(rng, 2, -5.0, 5.0);
        RowVec Ey = randomRow(rng, 2, -5.0, 5.0);
        RowVec Ez = randomRow(rng, 2, -5.0, 5.0);

        // Orientation vector perpendicular to the beam axis
        Eigen::Vector3d axis(Ex(1) - Ex(0), Ey(1) - Ey(0), Ez(1) - Ez(0));
        Eigen::Vector3d trial = randomRow(rng, 3, -1.0, 1.0).transpose();
        Eigen::Vector3d orient = axis.cross(trial);
        RowVec Eo = orient.transpose();

        RowVec Ep(6);
        Ep << 210e9, 80e9, 1e-2, 2e-5, 3e-5, 1e-5;
        RowVec RowEq = randomRow(rng, 4, -1e3, 1e3);

        Matrix Ke(12, 12), Ke_ref(12, 12);
        ColVec fe(12), fe_ref(12);

        beam3e(Ex, Ey, Ez, Eo, Ep, RowEq, Ke, fe);
        calfem_ref::beam3e(Ex, Ey, Ez, Eo, Ep, RowEq, Ke_ref, fe_ref);

        CHECK(Ke == Ke_ref);
        CHECK(fe == fe_ref);

        RowVec Ed = randomRow(rng, 12, -1e-3, 1e-3);
        int n = 2 + test % 19;

        Matrix Es, Edi, Es_ref, Edi_ref;
        ColVec Eci, Eci_ref;

        beam3s(Ex, Ey, Ez, Eo, Ep, Ed, RowEq, n, Es, Edi, Eci);
        calfem_ref::beam3s(Ex, Ey, Ez, Eo, Ep, Ed, RowEq, n, Es_ref, Edi_ref, Eci_ref);

        CHECK(Es == Es_ref);
        CHECK(Edi == Edi_ref);
        CHECK(Eci == Eci_ref);
    }
}

// ---- Full ObjectiveFrame-style solve against the reference ----

namespace {

// Space frame: four columns and four top beams (beam3e) plus two diagonal braces (bar3e),
// assembled the way ObjectiveFrame's BeamSolver does it.

struct SpaceFrame {
    static constexpr int nNodes = 8;
    Eigen::Index nDof = 6 * nNodes;

    Matrix coords;
    SpMatrix Ks;
    ColVec f;

    std::vector<std::array<int, 2>> beams{{0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
    std::vector<std::array<int, 2>> bars{{0, 5}, {2, 7}};

    RowVec Ep;
    RowVec Ep_b;

    SpaceFrame() : coords(nNodes, 3), Ep(6), Ep_b(2)
    {
        coords << 0, 0, 0, 4, 0, 0, 4, 3, 0, 0, 3, 0, //
            0, 0, 3, 4, 0, 3, 4, 3, 3, 0, 3, 3;

        Ep << 210e9, 80e9, 1e-2, 2e-5, 3e-5, 1e-5;
        Ep_b << 210e9, 5e-4;

        RowVec Ex(2), Ey(2), Ez(2), Eo(3);
        Matrix Ke(12, 12), Ke_b(6, 6);
        ColVec fe(12), fe_b(6);
        IntRowVec DofTopo(12), DofTopo_b(6);

        f.resize(nDof, 1);
        f.setZero();

        calfem::TripletList Ktriplets;

        for (auto &beam : beams)
        {
            elementCoords(beam, Ex, Ey, Ez);
            Eo = orientation(beam);

            RowVec RowEq(4);
            RowEq << 100.0, -2000.0, 500.0, 0.0;

            beam3e(Ex, Ey, Ez, Eo, Ep, RowEq, Ke, fe);
            beamTopo(beam, DofTopo);
            spassem(DofTopo, Ktriplets, Ke, f, fe);
        }

        for (auto &bar : bars)
        {
            elementCoords(bar, Ex, Ey, Ez);
            bar3e(Ex, Ey, Ez, Ep_b, -300.0, Ke_b, fe_b);
            barTopo(bar, DofTopo_b);
            spassem(DofTopo_b, Ktriplets, Ke_b, f, fe_b);
        }

        // Horizontal point loads at the top
        f(6 * 4 + 0) += 5e3;
        f(6 * 6 + 1) += -3e3;

        Ks.resize(nDof, nDof);
        Ks.setFromTriplets(Ktriplets.begin(), Ktriplets.end());
    }

    void elementCoords(const std::array<int, 2> &el, RowVec &Ex, RowVec &Ey, RowVec &Ez) const
    {
        Ex << coords(el[0], 0), coords(el[1], 0);
        Ey << coords(el[0], 1), coords(el[1], 1);
        Ez << coords(el[0], 2), coords(el[1], 2);
    }

    RowVec orientation(const std::array<int, 2> &el) const
    {
        RowVec Eo(3);
        if (coords(el[0], 2) != coords(el[1], 2))
            Eo << 1.0, 0.0, 0.0; // column
        else
            Eo << 0.0, 0.0, 1.0; // horizontal beam
        return Eo;
    }

    static void beamTopo(const std::array<int, 2> &el, IntRowVec &DofTopo)
    {
        for (int i = 0; i < 6; i++)
        {
            DofTopo(i) = 6 * el[0] + i + 1;
            DofTopo(6 + i) = 6 * el[1] + i + 1;
        }
    }

    static void barTopo(const std::array<int, 2> &el, IntRowVec &DofTopo_b)
    {
        for (int i = 0; i < 3; i++)
        {
            DofTopo_b(i) = 6 * el[0] + i + 1;
            DofTopo_b(3 + i) = 6 * el[1] + i + 1;
        }
    }

    // Fixed base, with a prescribed settlement of node 2 (sorted, as ObjectiveFrame builds them
    // from a std::map)
    static void baseSupports(IntColVec &bcDofs, ColVec &bcVals, double settlement)
    {
        bcDofs.resize(24);
        bcVals.resize(24);
        for (int i = 0; i < 24; i++)
        {
            bcDofs(i) = i + 1;
            bcVals(i) = 0.0;
        }
        bcVals(6 + 2) = settlement;
    }
};

} // namespace

TEST_CASE("ObjectiveFrame: SparseSolver setup/solve/recompute match reference", "[compat]")
{
    SpaceFrame model;

    std::vector<std::string> messages;
    setLogFunction([&messages](const std::string &message) { messages.push_back(message); });

    IntColVec m_bcDofs;
    ColVec m_bcVals;
    SpaceFrame::baseSupports(m_bcDofs, m_bcVals, -1e-3);

    // Solver objects held by value and reused, as in BeamSolver
    calfem::SparseSolver m_sparseSolver;
    calfem_ref::SparseSolver refSolver;

    ColVec m_globalA, m_globalQ, refA, refQ;

    m_globalA.resize(model.nDof, 1);
    m_globalA.setZero();
    refA = m_globalA;

    REQUIRE(m_sparseSolver.setup(model.Ks, m_bcDofs, m_bcVals));
    REQUIRE(refSolver.setup(model.Ks, m_bcDofs, m_bcVals));

    REQUIRE(m_sparseSolver.solve(model.f, m_globalA, m_globalQ));
    REQUIRE(refSolver.solve(model.f, refA, refQ));

    CHECK(m_globalA.isApprox(refA, 1e-12));
    CHECK(m_globalQ.isApprox(refQ, 1e-10));

    // Section forces from the solution, as in BeamSolver::execute
    for (auto &beam : model.beams)
    {
        RowVec Ex(2), Ey(2), Ez(2), Ed(12);
        model.elementCoords(beam, Ex, Ey, Ez);
        RowVec Eo = model.orientation(beam);

        IntRowVec DofTopo(12);
        SpaceFrame::beamTopo(beam, DofTopo);
        for (int i = 0; i < 12; i++)
            Ed(i) = m_globalA(DofTopo(i) - 1);

        RowVec RowEq(4);
        RowEq << 100.0, -2000.0, 500.0, 0.0;

        Matrix Es, Edi, Es_ref, Edi_ref;
        ColVec Eci, Eci_ref;
        beam3s(Ex, Ey, Ez, Eo, model.Ep, Ed, RowEq, 11, Es, Edi, Eci);
        calfem_ref::beam3s(Ex, Ey, Ez, Eo, model.Ep, Ed, RowEq, 11, Es_ref, Edi_ref, Eci_ref);

        CHECK(Es == Es_ref);
        CHECK(Edi == Edi_ref);
    }

    // Feedback-force recompute with a modified load vector
    ColVec f2 = model.f;
    f2(6 * 5 + 2) += -8e3;

    m_globalA.resize(model.nDof, 1);
    m_globalA.setZero();
    refA = m_globalA;

    REQUIRE(m_sparseSolver.recompute(f2, m_globalA, m_globalQ));
    REQUIRE(refSolver.recompute(f2, refA, refQ));

    CHECK(m_globalA.isApprox(refA, 1e-12));
    CHECK(m_globalQ.isApprox(refQ, 1e-10));

    // Re-running setup on the same solver objects with new boundary conditions
    SpaceFrame::baseSupports(m_bcDofs, m_bcVals, 2e-3);

    REQUIRE(m_sparseSolver.setup(model.Ks, m_bcDofs, m_bcVals));
    REQUIRE(refSolver.setup(model.Ks, m_bcDofs, m_bcVals));

    REQUIRE(m_sparseSolver.solve(model.f, m_globalA, m_globalQ));
    REQUIRE(refSolver.solve(model.f, refA, refQ));

    CHECK(m_globalA.isApprox(refA, 1e-12));
    CHECK(m_globalQ.isApprox(refQ, 1e-10));

    // A properly supported frame must not trigger the singularity warning
    setLogFunction(nullptr);
    for (auto &message : messages)
        CHECK(message.find("nearly singular") == std::string::npos);
}

TEST_CASE("ObjectiveFrame: spsolveq matches reference", "[compat]")
{
    SpaceFrame model;

    IntColVec bcDofs;
    ColVec bcVals;
    SpaceFrame::baseSupports(bcDofs, bcVals, -1e-3);

    ColVec a = ColVec::Zero(model.nDof), Q;
    ColVec refA = ColVec::Zero(model.nDof), refQ;

    REQUIRE(spsolveq(model.Ks, model.f, bcDofs, bcVals, a, Q));
    REQUIRE(calfem_ref::spsolveq(model.Ks, model.f, bcDofs, bcVals, refA, refQ));

    CHECK(a.isApprox(refA, 1e-12));
    CHECK(Q.isApprox(refQ, 1e-10));
}
