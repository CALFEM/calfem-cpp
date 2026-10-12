//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/assembly.h>

namespace calfem {

void assem(const IntRowVec &Topo, Matrix &K, const Matrix &Ke, ColVec &f, const ColVec &fe)
{
    int i, j, r, c;

    for (i = 0; i < Ke.rows(); i++)
        for (j = 0; j < Ke.cols(); j++)
        {
            r = Topo(i) - 1;
            c = Topo(j) - 1;

            K(r, c) = K(r, c) + Ke(i, j);
        }

    for (i = 0; i < fe.rows(); i++)
    {
        r = Topo(i) - 1;
        f(r) = f(r) + fe(i);
    }
}

void spassem(const IntRowVec &Topo, TripletList &Ktriplets, const Matrix &Ke, ColVec &f, const ColVec &fe)
{
    int i, j, r, c;

    for (i = 0; i < Ke.rows(); i++)
        for (j = 0; j < Ke.cols(); j++)
        {
            r = Topo(i) - 1;
            c = Topo(j) - 1;

            Ktriplets.push_back(Triplet(r, c, Ke(i, j)));
        }

    for (i = 0; i < fe.rows(); i++)
    {
        r = Topo(i) - 1;
        f(r) = f(r) + fe(i);
    }
}

// ---- Stiffness only ----

void assem(const IntRowVec &Topo, Matrix &K, const Matrix &Ke)
{
    for (Eigen::Index i = 0; i < Ke.rows(); i++)
        for (Eigen::Index j = 0; j < Ke.cols(); j++)
            K(Topo(i) - 1, Topo(j) - 1) += Ke(i, j);
}

void spassem(const IntRowVec &Topo, TripletList &Ktriplets, const Matrix &Ke)
{
    for (Eigen::Index i = 0; i < Ke.rows(); i++)
        for (Eigen::Index j = 0; j < Ke.cols(); j++)
            Ktriplets.push_back(Triplet(Topo(i) - 1, Topo(j) - 1, Ke(i, j)));
}

} // namespace calfem
