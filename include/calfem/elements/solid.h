//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Three-dimensional solid elements.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// ---- 8-node isoparametric brick, ep = [ir] (ir x ir x ir Gauss points, ir = 1, 2, 3) ----
//
// Three displacement dofs per node (ux, uy, uz). D is the 6 x 6 material matrix (hooke with
// ptype 4), eq = [bx, by, bz] a body force per unit volume (zero when omitted).
//
//   Ke = soli8e(ex, ey, ez, ep, D)            Matrix Ke = soli8e(ex, ey, ez, ep, D);
//   Ke, fe = soli8e(ex, ey, ez, ep, D, eq)    auto [Ke, fe] = soli8e(ex, ey, ez, ep, D, eq);
//   et, es, eci = soli8s(ex, ey, ez, ep, D, ed)
//                                             auto [et, es, eci] = soli8s(ex, ey, ez, ep, D, ed);
//
// NOTE: soli8s returns the STRAINS FIRST (et, es, eci), as calfem-python does. CALFEM for
// MATLAB returns [es, et, eci]. et and es hold [xx, yy, zz, xy, xz, yz] at each Gauss point
// (engineering shear strains), eci the Gauss point coordinates.

Matrix soli8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> soli8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                  const Matrix &D, const RowVec &eq);

std::tuple<Matrix, Matrix, Matrix> soli8s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                          const Matrix &D, const RowVec &ed);

} // namespace calfem
