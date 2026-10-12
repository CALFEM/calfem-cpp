//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Solvers for static FE equations with boundary conditions.

#pragma once

#include <calfem/types.h>

#include <Eigen/SparseCholesky>

#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace calfem {

// ---- Errors ----
//
// The value-returning solver functions (calfem-python style) throw SolverError when the
// input is inconsistent or the system cannot be solved. The output-argument versions
// return false instead. The message is the same one passed to the log function.

class SolverError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// ---- Solvers ----
//
// Solve K a = f with a(bcDofs) = bcVals. bcDofs are 1-based and may be given in any order;
// bcVals(i) is the value of dof bcDofs(i). a and Q (reaction forces, Q = K a - f) are resized.
// Returns false if the input is inconsistent or the factorization fails (K not positive definite
// on the free dofs, e.g. an unsupported structure).

bool solveq(const Matrix &K, const ColVec &f, const IntColVec &bcDofs, const ColVec &bcVals, ColVec &a, ColVec &Q);
bool spsolveq(const SpMatrix &K, const ColVec &f, const IntColVec &bcDofs, const ColVec &bcVals, ColVec &a, ColVec &Q);

// Value-returning versions, following calfem-python (a, r = solveq(K, f, bcPrescr, bcVal)).
// Without bcVals the prescribed values are zero. Throw SolverError on failure.
//
//   auto [a, r] = solveq(K, f, bcDofs, bcVals);
//   auto [a, r] = spsolveq(K, f, bcDofs);

std::tuple<ColVec, ColVec> solveq(const Matrix &K, const ColVec &f, const IntColVec &bcDofs);
std::tuple<ColVec, ColVec> solveq(const Matrix &K, const ColVec &f, const IntColVec &bcDofs, const ColVec &bcVals);
std::tuple<ColVec, ColVec> spsolveq(const SpMatrix &K, const ColVec &f, const IntColVec &bcDofs);
std::tuple<ColVec, ColVec> spsolveq(const SpMatrix &K, const ColVec &f, const IntColVec &bcDofs,
                                    const ColVec &bcVals);

// ---- Reduction, static condensation and eigenvalues (calfem-python style) ----
//
//   B = red(A, b)                 Matrix B = red(A, bcDofs);
//   K1, f1 = statcon(K, f, cd)    auto [K1, f1] = statcon(K, f, cd);
//   L, X = eigen(K, M[, b])       auto [L, X] = eigen(K, M, bcDofs);
//
// red removes the rows and columns of the (1-based) dofs b from A. b may also be a boundary
// condition matrix, in which case its first column holds the dofs.
//
// statcon condenses out the dofs cd: K1 = Kaa - Kab Kbb^-1 Kab^T, f1 = fa - Kab Kbb^-1 fb,
// where a are the remaining dofs (ascending) and b the dofs in cd. K must be symmetric.
//
// eigen solves K X = M X L for the eigenvalues L (ascending) and the mass-normalized
// eigenvectors X (one per column), with the dofs b (1-based) fixed. With b, X has ndof rows
// with zeros at the fixed dofs and one column per free dof. K must be symmetric and M
// symmetric positive definite (the signs of the eigenvectors are arbitrary). Throws
// SolverError if the problem cannot be solved.

Matrix red(const Matrix &A, const IntColVec &b);
Matrix red(const Matrix &A, const Matrix &bc);

std::tuple<Matrix, ColVec> statcon(const Matrix &K, const ColVec &f, const IntColVec &cd);

std::tuple<ColVec, Matrix> eigen(const Matrix &K, const Matrix &M);
std::tuple<ColVec, Matrix> eigen(const Matrix &K, const Matrix &M, const IntColVec &b);

// Sparse solver that keeps the reduced system and its factorization between solves.
//
//   setup()     partitions K into free/prescribed parts (copies what it needs; K, bcDofs and
//               bcVals do not have to outlive the solver).
//   solve()     factorizes the reduced system and solves for the load vector f.
//   recompute() solves for a new load vector reusing the factorization from solve()
//               (factorizes first if needed).
//
// solve() and recompute() also come in value-returning versions that throw SolverError:
//
//   auto [a, r] = solver.solve(f);

class SparseSolver {
private:
    Eigen::SimplicialLLT<SpMatrix> m_solver;
    bool m_isSetup{false};
    bool m_factorized{false};

    SpMatrix m_Ksys;  // free x free
    SpMatrix m_Ksysf; // free x prescribed
    SpMatrix m_Kp;    // prescribed rows of K (prescribed x all), for reaction forces

    Eigen::Index m_nFreeDofs{0};
    Eigen::Index m_nDofs{0};

    std::vector<Eigen::Index> m_ind;  // free dofs, zero-based, ascending
    std::vector<Eigen::Index> m_pind; // prescribed dofs, zero-based, in bcDofs order

    ColVec m_bcVals;

    bool factorize();
    bool solveSystem(const ColVec &f, ColVec &a, ColVec &Q);

public:
    SparseSolver();

    static std::shared_ptr<SparseSolver> create();

    bool setup(const SpMatrix &K, const IntColVec &bcDofs, const ColVec &bcVals);
    bool solve(const ColVec &f, ColVec &a, ColVec &Q);
    bool recompute(const ColVec &f, ColVec &a, ColVec &Q);

    bool setup(const SpMatrix &K, const IntColVec &bcDofs); // prescribed values zero
    std::tuple<ColVec, ColVec> solve(const ColVec &f);
    std::tuple<ColVec, ColVec> recompute(const ColVec &f);

    bool isFactorized() const;
    Eigen::Index nFreeDofs() const;
};

typedef std::shared_ptr<SparseSolver> SparseSolverPtr;

} // namespace calfem
