//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Plane stress and plane strain elements.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// Two displacement dofs per node (ux, uy). ptype = 1 plane stress, 2 plane strain. D is the
// material matrix from hooke: 3 x 3 ([xx, yy, xy]) or 4 x 4 ([xx, yy, zz, xy]); a 4 x 4 D is
// reduced to the plane case (for plane stress by condensing out the zz stress). eq = [bx, by]
// is a body force per unit volume (zero when omitted). Stresses es and strains et have as
// many components as D has columns.
//
//   Ke = plante(ex, ey, ep, D)              Matrix Ke = plante(ex, ey, ep, D);
//   Ke, fe = plante(ex, ey, ep, D, eq)      auto [Ke, fe] = plante(ex, ey, ep, D, eq);
//   es, et = plants(ex, ey, ep, D, ed)      auto [es, et] = plants(ex, ey, ep, D, ed);

// ---- 3-node triangle (constant strain), ep = [ptype, t] ----

Matrix plante(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> plante(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &eq);

// Stresses and strains, one row per row of ed. ex, ey hold one row per element, or a single
// row used for all rows of ed (several load cases of one element).
//
// Deviation from calfem-python: a single element with several rows in ed gives one result per
// row (calfem-python only evaluates the first row).
std::tuple<Matrix, Matrix> plants(const Matrix &ex, const Matrix &ey, const RowVec &ep, const Matrix &D,
                                  const Matrix &ed);

// Internal force vector from the element stresses es (3 or 4 components)
ColVec plantf(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &es);

// ---- 4-node quadrilateral built from four triangles, ep = [ptype, t] ----

Matrix planqe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> planqe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &eq);

// Area-weighted average of the stresses and strains of the four triangles (one element)
std::tuple<RowVec, RowVec> planqs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &ed);
std::tuple<RowVec, RowVec> planqs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const RowVec &ed, const RowVec &eq);

// ---- 4-node isoparametric quadrilateral, ep = [ptype, t, ir] (ir x ir Gauss points) ----
//
// Deviation from calfem-python: without eq only Ke is returned (calfem-python always
// returns Ke, fe).

Matrix plani4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> plani4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                   const RowVec &eq);

// Stresses and strains at the Gauss points (one row per point)
std::tuple<Matrix, Matrix> plani4s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                   const RowVec &ed);

} // namespace calfem
