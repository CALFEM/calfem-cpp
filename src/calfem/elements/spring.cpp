//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/elements/spring.h>

namespace calfem {

Matrix spring1e(double ep)
{
    double k = ep;

    Matrix Ke(2, 2);
    Ke << k, -k, //
        -k, k;
    return Ke;
}

double spring1s(double ep, const RowVec &ed)
{
    double k = ep;
    return k * (ed(1) - ed(0));
}

} // namespace calfem
