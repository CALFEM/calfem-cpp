//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Umbrella header: includes the complete CALFEM C++ library.

#pragma once

#include <calfem/types.h>
#include <calfem/log.h>
#include <calfem/utils.h>
#include <calfem/material.h>
#include <calfem/elements.h>
#include <calfem/dofs.h>
#include <calfem/assembly.h>
#include <calfem/boundary.h>
#include <calfem/solvers.h>
#include <calfem/dynamics.h>

#ifdef CALFEM_WITH_MESH
#include <calfem/geometry.h>
#include <calfem/mesh.h>
#endif

// Kept for source compatibility: code written against the original single-header library
// may rely on these being included.

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <functional>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
