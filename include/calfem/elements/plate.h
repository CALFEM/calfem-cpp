//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Plate elements.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// ---- Rectangular Kirchhoff plate element (4 nodes, 12 dofs) ----
//
//   Ke = platre(ex, ey, ep, D)          Matrix Ke = platre(ex, ey, ep, D);
//   Ke, fe = platre(ex, ey, ep, D, eq)  auto [Ke, fe] = platre(ex, ey, ep, D, eq);
//
// ex, ey: corner coordinates (counter-clockwise, sides parallel to the axes), ep = [t]
// (thickness), D: 3 x 3 plane stress material matrix (hooke with ptype 1), eq: distributed
// load per unit area. Dofs per node: [w, theta_x, theta_y].

Matrix platre(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> platre(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq);

} // namespace calfem
