//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Field problem (heat flow, groundwater flow, ...) elements.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// The elements solve -div(D grad(a)) = Q. D is the constitutive matrix (2 x 2 in 2D, 3 x 3 in
// 3D) and eq the distributed source Q per unit area/volume (zero when omitted). The *s
// routines take the nodal values ed with one row per load case (or element, flw2ts) and return
// the flux es = -D grad(a) and the gradient et = grad(a), one row per evaluation point.
//
//   Ke = flw2te(ex, ey, ep, D)              Matrix Ke = flw2te(ex, ey, ep, D);
//   Ke, fe = flw2te(ex, ey, ep, D, eq)      auto [Ke, fe] = flw2te(ex, ey, ep, D, eq);
//   es, et = flw2ts(ex, ey, D, ed)          auto [es, et] = flw2ts(ex, ey, D, ed);
//   es, et, eci = flw2i4s(ex, ey, ep, D, ed)
//
// Deviation from calfem-python: flw2te/flw2qe called with eq always return (Ke, fe), also
// for eq = 0 (calfem-python then returns only Ke).

// ---- 3-node triangle, ep = [t] (thickness) ----

Matrix flw2te(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> flw2te(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq);

// ex, ey, ed: one row per element; es, et: one row per element
std::tuple<Matrix, Matrix> flw2ts(const Matrix &ex, const Matrix &ey, const Matrix &D, const Matrix &ed);

// ---- 4-node quadrilateral built from four triangles, ep = [t] ----

Matrix flw2qe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> flw2qe(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq);

// Flux and gradient averaged over the four triangles; ed: one row per load case
std::tuple<Matrix, Matrix> flw2qs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const Matrix &ed);
std::tuple<Matrix, Matrix> flw2qs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                  const Matrix &ed, double eq);

// ---- 4-node isoparametric quadrilateral, ep = [t, ir] (ir x ir Gauss points, ir = 1, 2, 3) ----

Matrix flw2i4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> flw2i4e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq);

// es, et: one row per Gauss point and load case (Gauss point fastest); eci: Gauss point coordinates
std::tuple<Matrix, Matrix, Matrix> flw2i4s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                           const Matrix &ed);

// ---- 8-node isoparametric quadrilateral, ep = [t, ir] ----

Matrix flw2i8e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> flw2i8e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D, double eq);

std::tuple<Matrix, Matrix, Matrix> flw2i8s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &D,
                                           const Matrix &ed);

// ---- 8-node isoparametric brick, ep = [ir] (ir x ir x ir Gauss points, ir = 1, 2, 3) ----
//
// Deviation from calfem-python: ir = 3 uses the correct 27-point Gauss rule (in the CALFEM for
// MATLAB point order). calfem-python's rule for ir = 3 places all points in one zeta plane.

Matrix flw3i8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const Matrix &D);
std::tuple<Matrix, ColVec> flw3i8e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                   const Matrix &D, double eq);

std::tuple<Matrix, Matrix, Matrix> flw3i8s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                           const Matrix &D, const Matrix &ed);

} // namespace calfem
