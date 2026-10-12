//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Diagnostic logging.

#pragma once

#include <functional>
#include <string>

namespace calfem {

// ---- Logging ----
//
// Diagnostic messages from the solvers are passed to a user-supplied log function.
// No function is installed by default, so the library is silent.

using LogFunction = std::function<void(const std::string &message)>;

void setLogFunction(LogFunction func);

} // namespace calfem
