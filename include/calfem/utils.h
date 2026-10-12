//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Printing and inspection utilities.

#pragma once

#include <calfem/types.h>

#include <iostream>

namespace calfem {

template <typename T> void print(const T &A)
{
    std::cout << A.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]")) << "\n";
}

void spy(const SpMatrix &matrix);

// Scale factor for plotting displacements: rat times the largest element size divided by the
// largest displacement (ex, ey, ed: one row per element).
//
//   sfac = scalfact2(ex, ey, ed, rat)    double sfac = scalfact2(ex, ey, ed);
double scalfact2(const Matrix &ex, const Matrix &ey, const Matrix &ed, double rat = 0.2);

} // namespace calfem
