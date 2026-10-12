//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Degree-of-freedom handling (1-based dof numbers, as in CALFEM).

#pragma once

#include <calfem/types.h>

#include <tuple>
#include <vector>

namespace calfem {

std::vector<Eigen::Index> toZeroIndex(const std::vector<Eigen::Index> &vec);

// allDofs = 1..nDofs, freeDofs = allDofs without bcDofs (ascending). bcDofs may be unsorted.
void extractAllDofs(Eigen::Index nDofs, const std::vector<Eigen::Index> &bcDofs, std::vector<Eigen::Index> &freeDofs,
                    std::vector<Eigen::Index> &allDofs);

void extractElDisp(const IntMatrix &edof, const ColVec &a, Matrix &ed);

void extractEdof(const IntMatrix &dofs, const Matrix &coords, const IntMatrix &topo, IntMatrix &edof, Matrix &ex,
                 Matrix &ey, Matrix &ez);

// Value-returning versions:
//
//   auto [freeDofs, allDofs] = extractAllDofs(nDofs, bcDofs);
//   Matrix ed = extractElDisp(edof, a);                          (Python: ed = extract_eldisp(edof, a))
//   auto [edof, ex, ey, ez] = extractEdof(dofs, coords, topo);

std::tuple<std::vector<Eigen::Index>, std::vector<Eigen::Index>> extractAllDofs(
    Eigen::Index nDofs, const std::vector<Eigen::Index> &bcDofs);

Matrix extractElDisp(const IntMatrix &edof, const ColVec &a);

std::tuple<IntMatrix, Matrix, Matrix, Matrix> extractEdof(const IntMatrix &dofs, const Matrix &coords,
                                                          const IntMatrix &topo);

// ---- Dof numbering and element coordinates (calfem-python style) ----
//
//   dofs = create_dofs(nCoords, nDof)           IntMatrix dofs = create_dofs(nCoords, nDof);
//   ex, ey, ez = coord_extract(edof, coords, dofs[, nen])
//                                               auto [ex, ey, ez] = coord_extract(edof, coords, dofs);
//
// create_dofs numbers the dofs 1, 2, ... node by node (nCoords x nDof).
//
// coord_extract finds the coordinates of the element nodes from the element dofs (edof, one
// row per element), the node coordinates and the node dofs. nen is the number of nodes per
// element (default edof.cols() / dofs.cols()). Always returns ex, ey and ez; directions not
// present in coords are zero (calfem-python returns only as many as coords has columns).
// Throws std::invalid_argument if an element node is not found in dofs.

IntMatrix create_dofs(int nCoords, int nDof);

std::tuple<Matrix, Matrix, Matrix> coord_extract(const IntMatrix &edof, const Matrix &coords, const IntMatrix &dofs,
                                                 int nen = -1);

// Node coordinates, topology and node dofs from element data (for visualization):
//
//   coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node, ignore_first)
//
// Nodes are identified by their dofs and numbered in order of appearance. topo holds ZERO-BASED
// node numbers, as in calfem-python. With ignore_first, the first column of edof (an element
// number) is skipped.

std::tuple<Matrix, IntMatrix, IntMatrix> convert_to_node_topo(const IntMatrix &edof, const Matrix &ex,
                                                              const Matrix &ey, const Matrix &ez,
                                                              int n_dofs_per_node = 3, bool ignore_first = true);

// ---- Nodal values from element values ----
//
//   ev = stress2nodal(eseff, edof)    Matrix ev = stress2nodal(eseff, edof);
//
// eseff holds one value per element (e.g. effective stress), edof the element dofs with two dofs
// per node. Every node receives eseff / nnodes from each element it belongs to (as in
// calfem-python); ev holds the resulting value for each element node (one row per element).

Matrix stress2nodal(const ColVec &eseff, const IntMatrix &edof);

// calfem-python style names (aliases, both forms):
//
//   Matrix ed = extract_eldisp(edof, a);
//   auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);

inline Matrix extract_eldisp(const IntMatrix &edof, const ColVec &a)
{
    return extractElDisp(edof, a);
}

inline void extract_eldisp(const IntMatrix &edof, const ColVec &a, Matrix &ed)
{
    extractElDisp(edof, a, ed);
}

inline std::tuple<IntMatrix, Matrix, Matrix, Matrix> extract_edof(const IntMatrix &dofs, const Matrix &coords,
                                                                  const IntMatrix &topo)
{
    return extractEdof(dofs, coords, topo);
}

inline void extract_edof(const IntMatrix &dofs, const Matrix &coords, const IntMatrix &topo, IntMatrix &edof,
                         Matrix &ex, Matrix &ey, Matrix &ez)
{
    extractEdof(dofs, coords, topo, edof, ex, ey, ez);
}

} // namespace calfem
