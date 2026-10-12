//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Boundary conditions and loads on marked boundaries (calfem-python utils).

#pragma once

#include <calfem/types.h>

#include <map>
#include <tuple>
#include <vector>

namespace calfem {

// Dofs (1-based) on each marked boundary, marker -> dofs. With two dofs per node (2D) the dofs
// are stored node by node ([ux1, uy1, ux2, uy2, ...]), likewise with three dofs per node in 3D.
using BoundaryDofs = std::map<int, std::vector<int>>;

// ---- Prescribed values ----
//
//   bc, bcVal = apply_bc(bdofs, bc, bcVal, marker, value, dimension)
//   auto [bc, bcVal] = apply_bc(bdofs, bc, bcVal, marker, value, dimension);
//
// Adds the dofs of the boundary marker with the given value. dimension = 0 prescribes all dofs,
// 1 or 2 (1, 2 or 3 for apply_bc_3d) only the x, y (or z) dofs. The result is sorted and
// without duplicates; a dof that is already prescribed keeps its value. Throws
// std::invalid_argument for an unknown marker or dimension.

std::tuple<IntColVec, ColVec> apply_bc(const BoundaryDofs &boundaryDofs, const IntColVec &bcPrescr,
                                       const ColVec &bcVal, int marker, double value = 0.0, int dimension = 0);

std::tuple<IntColVec, ColVec> apply_bc_3d(const BoundaryDofs &boundaryDofs, const IntColVec &bcPrescr,
                                          const ColVec &bcVal, int marker, double value = 0.0, int dimension = 0);

// Prescribes the dofs of one node (nodeIdx is a zero-based row in dofs, as in calfem-python).
// dimension = 0: all dofs of the node, otherwise only dof number dimension. The dofs are
// appended (not sorted).
std::tuple<IntColVec, ColVec> apply_bc_node(int nodeIdx, const IntMatrix &dofs, const IntColVec &bcPrescr,
                                            const ColVec &bcVal, double value = 0.0, int dimension = 0);

// ---- Loads ----
//
// Add value to every selected dof of the boundary marker (apply_force, apply_force_3d) or
// distribute value evenly over them (apply_force_total, apply_force_total_3d). Throw
// std::invalid_argument for an unknown marker or dimension.

void apply_force(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value = 0.0, int dimension = 0);
void apply_force_3d(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value = 0.0, int dimension = 0);
void apply_force_total(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value = 0.0,
                       int dimension = 0);
void apply_force_total_3d(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value = 0.0,
                          int dimension = 0);

// Adds value to the dofs of one node (nodeIdx zero-based, dimension as in apply_bc_node).
//
// Deviation from calfem-python: the load is applied to the dof itself. calfem-python indexes f
// with the 1-based dof number and loads the following dof instead.
void apply_force_node(int nodeIdx, const IntMatrix &dofs, ColVec &f, double value = 0.0, int dimension = 0);

} // namespace calfem
