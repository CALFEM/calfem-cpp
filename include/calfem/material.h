//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Material models and stress measures.

#pragma once

#include <calfem/types.h>

namespace calfem {

// ---- Linear elastic isotropic material ----
//
//   D = hooke(ptype, E, v)        Matrix D = hooke(ptype, E, v);
//
// ptype 1: plane stress (3 x 3, components [xx, yy, xy])
//       2: plane strain (4 x 4, [xx, yy, zz, xy])
//       3: axisymmetry (4 x 4, [rr, zz, tt, rz])
//       4: three dimensional (6 x 6, [xx, yy, zz, xy, xz, yz])
//
// Throws std::invalid_argument for other values of ptype.

Matrix hooke(int ptype, double E, double v);

// ---- Effective (von Mises) stress ----
//
//   eseff = effmises(es, ptype)   ColVec eseff = effmises(es, ptype);
//
// es holds one stress state per row: ptype 1 (plane stress) 3 or 4 columns, ptype 2 and 3
// (plane strain, axisymmetry) 4 columns, ptype 4 (3D) 6 columns, ordered as in hooke.
// Throws std::invalid_argument if ptype or the number of columns does not match.

ColVec effmises(const Matrix &es, int ptype);

} // namespace calfem
