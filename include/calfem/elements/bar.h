//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Bar (truss) elements.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// ---- 3D bar element ----
//
// Output arguments are resized by the routines; pre-sizing them is allowed but not required.
//
// Each routine also has a calfem-python style form where results are returned:
//
//   Ke = bar3e(ex, ey, ez, ep)                 Matrix Ke = bar3e(ex, ey, ez, ep);
//   Ke, fe = bar3e(ex, ey, ez, ep, eq)         auto [Ke, fe] = bar3e(ex, ey, ez, ep, eq);
//   es = bar3s(ex, ey, ez, ep, ed)             ColVec es = bar3s(ex, ey, ez, ep, ed);
//   es, edi, eci = bar3s(..., ed, eq, nep)     auto [es, edi, eci] = bar3s(..., ed, eq, nep);
//
// Without eq the load is zero. Without nep the section forces are evaluated at the two
// element ends, and only es is returned.

void bar3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, double eq, Matrix &Ke, ColVec &fe);

void bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const RowVec &ed, double eq, int n,
           ColVec &es, ColVec &edi, ColVec &eci);

Matrix bar3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep);
std::tuple<Matrix, ColVec> bar3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, double eq);

ColVec bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const RowVec &ed);
ColVec bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep, const RowVec &ed, double eq);
std::tuple<ColVec, ColVec, ColVec> bar3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &ep,
                                         const RowVec &ed, double eq, int nep);

// ---- 1D and 2D bar elements (calfem-python style) ----
//
// eq is the distributed axial load qX (zero when omitted). The *s routines evaluate the
// section forces at nep equally spaced points along the element (default 2: the element
// ends); without nep only es is returned. nep must be at least 2.
//
//   Ke = bar1e(ex, ep)                         Matrix Ke = bar1e(ex, ep);
//   Ke, fe = bar1e(ex, ep, eq)                 auto [Ke, fe] = bar1e(ex, ep, eq);
//   es = bar1s(ex, ep, ed)                     ColVec es = bar1s(ex, ep, ed);
//   es, edi, eci = bar1s(ex, ep, ed, eq, nep)  auto [es, edi, eci] = bar1s(ex, ep, ed, eq, nep);
//
// bar1e/bar1s:   ex = [x1, x2], ep = [E, A], ed = [u1, u2]
// bar1we/bar1ws: as bar1e/bar1s on an elastic foundation, ep = [E, A, kX]
// bar2e/bar2s:   ex = [x1, x2], ey = [y1, y2], ep = [E, A], ed = [u1, u2, u3, u4]
// bar2ge:        geometric nonlinear (second order) bar with axial force QX
// bar2gs:        returns es and the axial force QX (as Python)

Matrix bar1e(const RowVec &ex, const RowVec &ep);
std::tuple<Matrix, ColVec> bar1e(const RowVec &ex, const RowVec &ep, double eq);

ColVec bar1s(const RowVec &ex, const RowVec &ep, const RowVec &ed);
ColVec bar1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq);
std::tuple<ColVec, ColVec, ColVec> bar1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep);

Matrix bar1we(const RowVec &ex, const RowVec &ep);
std::tuple<Matrix, ColVec> bar1we(const RowVec &ex, const RowVec &ep, double eq);

ColVec bar1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed);
ColVec bar1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq);
std::tuple<ColVec, ColVec, ColVec> bar1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep);

Matrix bar2e(const RowVec &ex, const RowVec &ey, const RowVec &ep);
std::tuple<Matrix, ColVec> bar2e(const RowVec &ex, const RowVec &ey, const RowVec &ep, double eq);

ColVec bar2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed);
ColVec bar2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, double eq);
std::tuple<ColVec, ColVec, ColVec> bar2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                         double eq, int nep);

Matrix bar2ge(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX);

std::tuple<ColVec, double> bar2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed);
std::tuple<ColVec, double, ColVec, ColVec> bar2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep,
                                                  const RowVec &ed, int nep);

} // namespace calfem
