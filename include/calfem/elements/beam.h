//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Beam elements.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// ---- 3D beam element ----
//
// Output arguments are resized by the routines; pre-sizing them is allowed but not required.
//
// Each routine also has a calfem-python style form where results are returned:
//
//   Ke = beam3e(ex, ey, ez, eo, ep)            Matrix Ke = beam3e(ex, ey, ez, eo, ep);
//   Ke, fe = beam3e(ex, ey, ez, eo, ep, eq)    auto [Ke, fe] = beam3e(ex, ey, ez, eo, ep, eq);
//   es = beam3s(..., ed)                       Matrix es = beam3s(..., ed);
//   es, edi, eci = beam3s(..., ed, eq, nep)    auto [es, edi, eci] = beam3s(..., ed, eq, nep);
//
// Without eq the load is zero. Without nep the section forces are evaluated at the two
// element ends, and only es is returned.

void beam3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep, const RowVec &eq,
            Matrix &Ke, ColVec &fe);

void beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep, const RowVec &ed,
            const RowVec &eq, const int n, Matrix &es, Matrix &edi, ColVec &eci);

Matrix beam3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep);
std::tuple<Matrix, ColVec> beam3e(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo,
                                  const RowVec &ep, const RowVec &eq);

Matrix beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep,
              const RowVec &ed);
Matrix beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo, const RowVec &ep,
              const RowVec &ed, const RowVec &eq);
std::tuple<Matrix, Matrix, ColVec> beam3s(const RowVec &ex, const RowVec &ey, const RowVec &ez, const RowVec &eo,
                                          const RowVec &ep, const RowVec &ed, const RowVec &eq, int nep);

// ---- 1D beam elements (calfem-python style) ----
//
// eq is the distributed transverse load qY (zero when omitted). The *s routines evaluate
// at nep equally spaced points (default 2, minimum 2); without nep only es is returned.
//
//   Ke = beam1e(ex, ep)                        Matrix Ke = beam1e(ex, ep);
//   Ke, fe = beam1e(ex, ep, eq)                auto [Ke, fe] = beam1e(ex, ep, eq);
//   es = beam1s(ex, ep, ed)                    Matrix es = beam1s(ex, ep, ed);
//   es, edi, eci = beam1s(ex, ep, ed, eq, nep) auto [es, edi, eci] = beam1s(ex, ep, ed, eq, nep);
//
// beam1e/beam1s:   ex = [x1, x2], ep = [E, I], ed = [v1, theta1, v2, theta2],
//                  es = [V, M] per point, edi = v
// beam1we/beam1ws: on an elastic foundation, ep = [E, I, kY]

Matrix beam1e(const RowVec &ex, const RowVec &ep);
std::tuple<Matrix, ColVec> beam1e(const RowVec &ex, const RowVec &ep, double eq);

Matrix beam1s(const RowVec &ex, const RowVec &ep, const RowVec &ed);
Matrix beam1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq);
std::tuple<Matrix, ColVec, ColVec> beam1s(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep);

Matrix beam1we(const RowVec &ex, const RowVec &ep);
std::tuple<Matrix, ColVec> beam1we(const RowVec &ex, const RowVec &ep, double eq);

Matrix beam1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed);
Matrix beam1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq);
std::tuple<Matrix, ColVec, ColVec> beam1ws(const RowVec &ex, const RowVec &ep, const RowVec &ed, double eq, int nep);

// ---- 2D beam elements (calfem-python style) ----
//
// ex = [x1, x2], ey = [y1, y2], ed = [u1, ..., u6] (two nodes with ux, uy, theta).
// eq = [qX, qY] (local axial and transverse distributed load, zero when omitted).
// es = [N, V, M] per evaluation point, edi = [u, v] (beam2ts: [u, v, theta]).
//
// beam2e/beam2s:     ep = [E, A, I]
// beam2we/beam2ws:   on an elastic foundation, ep = [E, A, I, kX, kY]
// beam2te/beam2ts:   Timoshenko beam, ep = [E, G, A, I, ks]
// beam2ge/beam2gs:   second order theory with axial force QX, ep = [E, A, I], eq = qY.
//                    beam2gs also returns the updated axial force QX.
// beam2gxe/beam2gxs: as beam2ge/beam2gs using the exact solution of the differential
//                    equation instead of polynomial shape functions.
// beam2de/beam2ds:   dynamic analysis, ep = [E, A, I, m] or [E, A, I, m, a, b] (Rayleigh
//                    damping C = a M + b K). beam2de always returns (Ke, Me, Ce); Ce is
//                    zero without damping (Python returns only Ke, Me in that case).
//                    beam2ds takes one row per time step in ed, ev and ea.

Matrix beam2e(const RowVec &ex, const RowVec &ey, const RowVec &ep);
std::tuple<Matrix, ColVec> beam2e(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &eq);

Matrix beam2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed);
Matrix beam2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, const RowVec &eq);
std::tuple<Matrix, Matrix, ColVec> beam2s(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                          const RowVec &eq, int nep);

Matrix beam2we(const RowVec &ex, const RowVec &ey, const RowVec &ep);
std::tuple<Matrix, ColVec> beam2we(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &eq);

Matrix beam2ws(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed);
Matrix beam2ws(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, const RowVec &eq);
std::tuple<Matrix, Matrix, ColVec> beam2ws(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                           const RowVec &eq, int nep);

Matrix beam2te(const RowVec &ex, const RowVec &ey, const RowVec &ep);
std::tuple<Matrix, ColVec> beam2te(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &eq);

Matrix beam2ts(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed);
Matrix beam2ts(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed, const RowVec &eq);
std::tuple<Matrix, Matrix, ColVec> beam2ts(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                           const RowVec &eq, int nep);

Matrix beam2ge(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX);
std::tuple<Matrix, ColVec> beam2ge(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX, double eq);

std::tuple<Matrix, double> beam2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                   double QX);
std::tuple<Matrix, double> beam2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                   double QX, double eq);
std::tuple<Matrix, double, Matrix, ColVec> beam2gs(const RowVec &ex, const RowVec &ey, const RowVec &ep,
                                                   const RowVec &ed, double QX, double eq, int nep);

Matrix beam2gxe(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX);
std::tuple<Matrix, ColVec> beam2gxe(const RowVec &ex, const RowVec &ey, const RowVec &ep, double QX, double eq);

std::tuple<Matrix, double> beam2gxs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                    double QX);
std::tuple<Matrix, double> beam2gxs(const RowVec &ex, const RowVec &ey, const RowVec &ep, const RowVec &ed,
                                    double QX, double eq);
std::tuple<Matrix, double, Matrix, ColVec> beam2gxs(const RowVec &ex, const RowVec &ey, const RowVec &ep,
                                                    const RowVec &ed, double QX, double eq, int nep);

std::tuple<Matrix, Matrix, Matrix> beam2de(const RowVec &ex, const RowVec &ey, const RowVec &ep);

Matrix beam2ds(const RowVec &ex, const RowVec &ey, const RowVec &ep, const Matrix &ed, const Matrix &ev,
               const Matrix &ea);

// ---- Deformed shape of 2D beams (for plotting) ----
//
//   excd, eycd = beam2crd(ex, ey, ed, mag)    auto [excd, eycd] = beam2crd(ex, ey, ed, mag);
//
// ex, ey, ed: one row per element. Returns the coordinates of 21 points along each deformed
// element (one row per element), with the displacements scaled by mag.

std::tuple<Matrix, Matrix> beam2crd(const Matrix &ex, const Matrix &ey, const Matrix &ed, double mag);

} // namespace calfem
