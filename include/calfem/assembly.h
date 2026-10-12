//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Assembly of element matrices into global matrices.

#pragma once

#include <calfem/types.h>

namespace calfem {

void assem(const IntRowVec &Topo, Matrix &K, const Matrix &Ke, ColVec &f, const ColVec &fe);

void spassem(const IntRowVec &Topo, TripletList &Ktriplets, const Matrix &Ke, ColVec &f, const ColVec &fe);

// Stiffness only, as calfem-python's assem(edof, K, Ke) without f and fe. K (or the triplet
// list) is updated in place, as returning it would copy the global matrix for every element.

void assem(const IntRowVec &Topo, Matrix &K, const Matrix &Ke);

void spassem(const IntRowVec &Topo, TripletList &Ktriplets, const Matrix &Ke);

} // namespace calfem
