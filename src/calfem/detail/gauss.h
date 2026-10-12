//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Gauss integration rules for isoparametric elements (not installed).

#pragma once

#include <calfem/types.h>

#include <stdexcept>
#include <string>

namespace calfem::detail {

struct GaussRule {
    Matrix points;  // ngp x dim, natural coordinates
    ColVec weights; // ngp, products of the 1D weights
};

// The CALFEM constants (as in CALFEM for MATLAB and calfem-python)
constexpr double gaussG2 = 0.577350269189626;
constexpr double gaussG3 = 0.774596669241483;
constexpr double gaussW3a = 0.555555555555555;
constexpr double gaussW3b = 0.888888888888888;

// 2D rule with ir x ir points, ir = 1, 2 or 3, in CALFEM order (xsi fastest).

inline GaussRule gaussRule2D(int ir, const char *caller)
{
    GaussRule r;

    if (ir == 1)
    {
        r.points = Matrix::Zero(1, 2);
        r.weights = ColVec::Constant(1, 4.0);
    }
    else if ((ir == 2) || (ir == 3))
    {
        double p2[] = {-gaussG2, gaussG2};
        double w2[] = {1.0, 1.0};
        double p3[] = {-gaussG3, 0.0, gaussG3};
        double w3[] = {gaussW3a, gaussW3b, gaussW3a};

        const double *p = ir == 2 ? p2 : p3;
        const double *w = ir == 2 ? w2 : w3;

        r.points.resize(ir * ir, 2);
        r.weights.resize(ir * ir);

        for (int j = 0; j < ir; j++)
            for (int i = 0; i < ir; i++)
            {
                int k = i + ir * j;
                r.points(k, 0) = p[i];
                r.points(k, 1) = p[j];
                r.weights(k) = w[i] * w[j];
            }
    }
    else
        throw std::invalid_argument(std::string(caller) + ": the integration rule ir must be 1, 2 or 3.");

    return r;
}

// 3D rule with ir x ir x ir points, ir = 1, 2 or 3.
//
// ir = 2 uses the CALFEM point order (counter-clockwise in each zeta plane), ir = 3 the
// CALFEM for MATLAB order (xsi fastest, then eta, then zeta).

inline GaussRule gaussRule3D(int ir, const char *caller)
{
    GaussRule r;

    if (ir == 1)
    {
        r.points = Matrix::Zero(1, 3);
        r.weights = ColVec::Constant(1, 8.0);
    }
    else if (ir == 2)
    {
        r.points.resize(8, 3);
        r.points << -1, -1, -1, //
            1, -1, -1,          //
            1, 1, -1,           //
            -1, 1, -1,          //
            -1, -1, 1,          //
            1, -1, 1,           //
            1, 1, 1,            //
            -1, 1, 1;
        r.points *= gaussG2;
        r.weights = ColVec::Ones(8);
    }
    else if (ir == 3)
    {
        double p[] = {-gaussG3, 0.0, gaussG3};
        double w[] = {gaussW3a, gaussW3b, gaussW3a};

        r.points.resize(27, 3);
        r.weights.resize(27);

        for (int k = 0; k < 3; k++)
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 3; i++)
                {
                    int n = i + 3 * j + 9 * k;
                    r.points(n, 0) = p[i];
                    r.points(n, 1) = p[j];
                    r.points(n, 2) = p[k];
                    r.weights(n) = w[i] * w[j] * w[k];
                }
    }
    else
        throw std::invalid_argument(std::string(caller) + ": the integration rule ir must be 1, 2 or 3.");

    return r;
}

} // namespace calfem::detail
