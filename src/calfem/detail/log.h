//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Internal logging helpers shared by the library sources (not installed).

#pragma once

#include <string>

namespace calfem::detail {

// Passes message to the user log function, if one is installed.
void logMessage(const std::string &message);

// Logs message and remembers it as the last error of the calling thread.
void logError(const std::string &message);

// Last error message of the calling thread (used for SolverError messages).
std::string &lastError();

} // namespace calfem::detail
