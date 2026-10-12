//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Spring elements.

#pragma once

#include <calfem/types.h>

namespace calfem {

// ---- 1D spring element ----
//
//   Ke = spring1e(ep)        Matrix Ke = spring1e(k);
//   es = spring1s(ep, ed)    double N = spring1s(k, ed);
//
// ep = k (spring stiffness), ed = [u1, u2].

Matrix spring1e(double ep);

double spring1s(double ep, const RowVec &ed);

} // namespace calfem
