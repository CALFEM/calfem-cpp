//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// CALFEM plotpar colour and line type codes, shared by the vis sources (not installed).

#pragma once

#include <calfem/vis.h>

#include <stdexcept>
#include <string>

namespace calfem::vis::codes {

// Colour code: 1 black, 2 blue, 3 magenta, 4 red.
inline Color colorCode(int code, const std::string &what)
{
    switch (code)
    {
    case 1:
        return {0.0, 0.0, 0.0};
    case 2:
        return {0.0, 0.0, 1.0};
    case 3:
        return {1.0, 0.0, 1.0};
    case 4:
        return {1.0, 0.0, 0.0};
    default:
        throw std::invalid_argument("Invalid value for " + what + " (colour code 1-4), got " + std::to_string(code) +
                                    ".");
    }
}

// Line type code: 1 solid, 2 dashed, 3 dotted.
inline LineStyle lineCode(int code, const std::string &what)
{
    switch (code)
    {
    case 1:
        return LineStyle::Solid;
    case 2:
        return LineStyle::Dashed;
    case 3:
        return LineStyle::Dotted;
    default:
        throw std::invalid_argument("Invalid value for " + what + " (line type 1-3), got " + std::to_string(code) + ".");
    }
}

} // namespace calfem::vis::codes
