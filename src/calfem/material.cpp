//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Material models, ported from calfem-python (core.py).

#include <calfem/material.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace calfem {

Matrix hooke(int ptype, double E, double v)
{
    if (ptype == 1)
    {
        Matrix D(3, 3);
        D << 1.0, v, 0.0, //
            v, 1.0, 0.0,  //
            0.0, 0.0, (1.0 - v) / 2.0;
        return E * D / (1.0 - v * v);
    }
    else if ((ptype == 2) || (ptype == 3))
    {
        Matrix D(4, 4);
        D << 1.0 - v, v, v, 0.0,  //
            v, 1.0 - v, v, 0.0,   //
            v, v, 1.0 - v, 0.0,   //
            0.0, 0.0, 0.0, (1.0 - 2.0 * v) / 2.0;
        return E / (1.0 + v) * D / (1.0 - 2.0 * v);
    }
    else if (ptype == 4)
    {
        Matrix D = Matrix::Zero(6, 6);
        D.topLeftCorner(3, 3) << 1.0 - v, v, v, //
            v, 1.0 - v, v,                     //
            v, v, 1.0 - v;
        D(3, 3) = D(4, 4) = D(5, 5) = (1.0 - 2.0 * v) / 2.0;
        return E * D / (1.0 + v) / (1.0 - 2.0 * v);
    }

    throw std::invalid_argument("hooke: ptype must be 1, 2, 3 or 4, got " + std::to_string(ptype) + ".");
}

ColVec effmises(const Matrix &es, int ptype)
{
    auto ncomps = es.cols();
    ColVec seff(es.rows());

    if (ptype == 1)
    {
        if ((ncomps != 3) && (ncomps != 4))
            throw std::invalid_argument("effmises: es must have 3 or 4 columns for plane stress (ptype 1).");

        auto ixy = ncomps == 3 ? 2 : 3;
        for (Eigen::Index i = 0; i < es.rows(); i++)
        {
            double sxx = es(i, 0), syy = es(i, 1), sxy = es(i, ixy);
            seff(i) = std::sqrt(sxx * sxx + syy * syy - sxx * syy + 3.0 * sxy * sxy);
        }
    }
    else if ((ptype == 2) || (ptype == 3))
    {
        if (ncomps != 4)
            throw std::invalid_argument("effmises: es must have 4 columns for plane strain and axisymmetry (ptype 2, 3).");

        for (Eigen::Index i = 0; i < es.rows(); i++)
        {
            double s1 = es(i, 0), s2 = es(i, 1), s3 = es(i, 2), t12 = es(i, 3);
            seff(i) = std::sqrt(s1 * s1 + s2 * s2 + s3 * s3 - s1 * s2 - s2 * s3 - s3 * s1 + 3.0 * t12 * t12);
        }
    }
    else if (ptype == 4)
    {
        if (ncomps != 6)
            throw std::invalid_argument("effmises: es must have 6 columns for three dimensional analysis (ptype 4).");

        for (Eigen::Index i = 0; i < es.rows(); i++)
        {
            double sx = es(i, 0), sy = es(i, 1), sz = es(i, 2);
            double txy = es(i, 3), txz = es(i, 4), tyz = es(i, 5);
            seff(i) = std::sqrt(sx * sx + sy * sy + sz * sz - sx * sy - sy * sz - sz * sx +
                                3.0 * (txy * txy + txz * txz + tyz * tyz));
        }
    }
    else
        throw std::invalid_argument("effmises: ptype must be 1, 2, 3 or 4, got " + std::to_string(ptype) + ".");

    return seff;
}

} // namespace calfem
