//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Shape functions for the isoparametric elements (not installed).

#pragma once

#include <calfem/types.h>

namespace calfem::detail {

// Isoparametric shape functions: N (1 x nnodes) and dN/d(xsi, eta, zeta) (dim x nnodes),
// node order as in CALFEM.

inline void quad4Shape(double xsi, double eta, RowVec &N, Matrix &dNr)
{
    N.resize(4);
    N << (1 - xsi) * (1 - eta) / 4.0, (1 + xsi) * (1 - eta) / 4.0, (1 + xsi) * (1 + eta) / 4.0,
        (1 - xsi) * (1 + eta) / 4.0;

    dNr.resize(2, 4);
    dNr << -(1 - eta) / 4.0, (1 - eta) / 4.0, (1 + eta) / 4.0, -(1 + eta) / 4.0, //
        -(1 - xsi) / 4.0, -(1 + xsi) / 4.0, (1 + xsi) / 4.0, (1 - xsi) / 4.0;
}

inline void quad8Shape(double xsi, double eta, RowVec &N, Matrix &dNr)
{
    N.resize(8);
    N << -(1 - xsi) * (1 - eta) * (1 + xsi + eta) / 4.0, //
        -(1 + xsi) * (1 - eta) * (1 - xsi + eta) / 4.0,  //
        -(1 + xsi) * (1 + eta) * (1 - xsi - eta) / 4.0,  //
        -(1 - xsi) * (1 + eta) * (1 + xsi - eta) / 4.0,  //
        (1 - xsi * xsi) * (1 - eta) / 2.0,               //
        (1 + xsi) * (1 - eta * eta) / 2.0,               //
        (1 - xsi * xsi) * (1 + eta) / 2.0,               //
        (1 - xsi) * (1 - eta * eta) / 2.0;

    dNr.resize(2, 8);
    dNr << -(-(1 - eta) * (1 + xsi + eta) + (1 - xsi) * (1 - eta)) / 4.0, //
        -((1 - eta) * (1 - xsi + eta) - (1 + xsi) * (1 - eta)) / 4.0,     //
        -((1 + eta) * (1 - xsi - eta) - (1 + xsi) * (1 + eta)) / 4.0,     //
        -(-(1 + eta) * (1 + xsi - eta) + (1 - xsi) * (1 + eta)) / 4.0,    //
        -xsi * (1 - eta),                                                 //
        (1 - eta * eta) / 2.0,                                            //
        -xsi * (1 + eta),                                                 //
        -(1 - eta * eta) / 2.0,                                           //
        -(-(1 - xsi) * (1 + xsi + eta) + (1 - xsi) * (1 - eta)) / 4.0,    //
        -(-(1 + xsi) * (1 - xsi + eta) + (1 + xsi) * (1 - eta)) / 4.0,    //
        -((1 + xsi) * (1 - xsi - eta) - (1 + xsi) * (1 + eta)) / 4.0,     //
        -((1 - xsi) * (1 + xsi - eta) - (1 - xsi) * (1 + eta)) / 4.0,     //
        -(1 - xsi * xsi) / 2.0,                                           //
        -eta * (1 + xsi),                                                 //
        (1 - xsi * xsi) / 2.0,                                            //
        -eta * (1 - xsi);
}

inline void brick8Shape(double xsi, double eta, double zet, RowVec &N, Matrix &dNr)
{
    N.resize(8);
    N << (1 - xsi) * (1 - eta) * (1 - zet) / 8.0, (1 + xsi) * (1 - eta) * (1 - zet) / 8.0,
        (1 + xsi) * (1 + eta) * (1 - zet) / 8.0, (1 - xsi) * (1 + eta) * (1 - zet) / 8.0,
        (1 - xsi) * (1 - eta) * (1 + zet) / 8.0, (1 + xsi) * (1 - eta) * (1 + zet) / 8.0,
        (1 + xsi) * (1 + eta) * (1 + zet) / 8.0, (1 - xsi) * (1 + eta) * (1 + zet) / 8.0;

    dNr.resize(3, 8);
    dNr << -(1 - eta) * (1 - zet), (1 - eta) * (1 - zet), (1 + eta) * (1 - zet), -(1 + eta) * (1 - zet),
        -(1 - eta) * (1 + zet), (1 - eta) * (1 + zet), (1 + eta) * (1 + zet), -(1 + eta) * (1 + zet), //
        -(1 - xsi) * (1 - zet), -(1 + xsi) * (1 - zet), (1 + xsi) * (1 - zet), (1 - xsi) * (1 - zet),
        -(1 - xsi) * (1 + zet), -(1 + xsi) * (1 + zet), (1 + xsi) * (1 + zet), (1 - xsi) * (1 + zet), //
        -(1 - xsi) * (1 - eta), -(1 + xsi) * (1 - eta), -(1 + xsi) * (1 + eta), -(1 - xsi) * (1 + eta),
        (1 - xsi) * (1 - eta), (1 + xsi) * (1 - eta), (1 + xsi) * (1 + eta), (1 - xsi) * (1 + eta);
    dNr /= 8.0;
}

} // namespace calfem::detail
