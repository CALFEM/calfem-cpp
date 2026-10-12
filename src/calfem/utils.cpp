//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/utils.h>

#include <cmath>
#include <iostream>

#include <algorithm>

namespace calfem {

void spy(const SpMatrix &matrix)
{
    auto rows = matrix.rows();
    auto cols = matrix.cols();

    for (auto i = 0; i < rows; i++)
    {
        for (auto j = 0; j < cols; j++)
        {
            if (std::abs(matrix.coeff(i, j)) > 0.0)
                std::cout << "o ";
            else
                std::cout << ". ";
        }
        std::cout << "\n";
    }
}

double scalfact2(const Matrix &ex, const Matrix &ey, const Matrix &ed, double rat)
{
    double dxmax = (ex.rowwise().maxCoeff() - ex.rowwise().minCoeff()).maxCoeff();
    double dymax = (ey.rowwise().maxCoeff() - ey.rowwise().minCoeff()).maxCoeff();
    double dlmax = std::max(dxmax, dymax);
    double edmax = ed.cwiseAbs().maxCoeff();

    return rat * dlmax / edmax;
}

} // namespace calfem
