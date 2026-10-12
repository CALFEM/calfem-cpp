//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Internal interface between the mesh front end and the mesh generators (not installed).

#pragma once

#include <calfem/mesh.h>

#include <Eigen/Dense>

#include <map>
#include <set>
#include <vector>

namespace calfem::detail {

// Mesh as produced by a generator, before dofs are assigned.
struct RawMesh {
    int dim = 2;
    std::vector<Eigen::Vector3d> nodes;              // node coordinates
    std::vector<std::vector<int>> elements;          // element nodes, 1-based
    std::vector<int> markers;                        // element markers
    std::map<int, std::set<int>> boundaryNodes;      // marker -> nodes (1-based)
    BoundaryElements boundaryElements;               // marker -> boundary elements
};

RawMesh meshBuiltIn(const Geometry &geometry, int elType, double elSizeFactor);

// Defined in gmsh_mesher.cpp; throws if calfem was built without Gmsh
RawMesh meshGmsh(const Geometry &geometry, const MeshGenerator &generator);

// Assigns dofs and builds the Mesh (edof, dofs, bdofs)
Mesh finalizeMesh(const RawMesh &raw, int dofsPerNode);

} // namespace calfem::detail
