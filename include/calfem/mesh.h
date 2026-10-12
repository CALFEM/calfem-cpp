//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Mesh generation from a Geometry (port of calfem-python's calfem.mesh).

#pragma once

#include <calfem/boundary.h>
#include <calfem/geometry.h>
#include <calfem/types.h>

#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace calfem {

// A boundary element (line or point element on a marked curve/point), node numbers 1-based.
// elm_type uses the Gmsh numbering (1: 2-node line, 15: point).
struct BoundaryElement {
    int elm_type = 1;
    std::vector<int> nodes;
};

using BoundaryElements = std::map<int, std::vector<BoundaryElement>>;

struct Mesh {
    Matrix coords;                       // node coordinates, one row per node (2 or 3 columns)
    IntMatrix edof;                      // element dofs, one row per element
    IntMatrix dofs;                      // dofs of each node (1-based)
    BoundaryDofs bdofs;                  // marker -> dofs on the marked curves/points (2D) or surfaces (3D)
    std::vector<int> element_markers;    // marker of the surface (volume) of each element
    BoundaryElements boundary_elements;  // marker -> boundary elements
    IntMatrix topo;                      // element node numbers (1-based), one row per element
};

// Which mesh generator to use.
//
//   BuiltIn  the library's own 2D mesher (structured surfaces, CDT based triangulation)
//   Gmsh     the Gmsh library (requires calfem built with CALFEM_WITH_GMSH)
//   Auto     Gmsh if available, otherwise BuiltIn
enum class MeshBackend { Auto, BuiltIn, Gmsh };

// True if the library was built with Gmsh support (CALFEM_WITH_GMSH).
bool gmsh_available();

// Mesh generator, as calfem-python's GmshMeshGenerator:
//
//   MeshGenerator mesh(g);
//   mesh.el_type = 3;
//   mesh.el_size_factor = 0.05;
//   mesh.dofs_per_node = 2;
//   Mesh m = mesh.create();
//
// el_type uses the Gmsh element numbering: 2 = 3-node triangles, 3 = 4-node quadrilaterals,
// 4 = 4-node tetrahedra, 5 = 8-node hexahedra (and higher order types with Gmsh). Element sizes
// are the point el_size values multiplied by el_size_factor.
//
// The built-in mesher handles 2D geometries with el_type 2 or 3:
//   - curves are divided according to el_on_curve or the element sizes at their end points
//   - structured surfaces (struct_surface) get a transfinite (mapped) mesh
//   - other surfaces are triangulated (constrained Delaunay, graded interior points, smoothing);
//     for el_type 3 every triangle is split into three quadrilaterals, which requires an even
//     number of elements on the curves of such surfaces
// 3D geometries, other element types and min_size/max_size/meshing_algorithm/clcurv require Gmsh.
//
// The Gmsh backend writes the geometry as a .geo file in the same way as calfem-python, meshes it
// with the Gmsh library and returns the same data as calfem-python.

class MeshGenerator {
public:
    explicit MeshGenerator(const Geometry &geometry);

    int el_type = 2;
    double el_size_factor = 1.0;
    int dofs_per_node = 1;

    MeshBackend backend = MeshBackend::Auto;

    // Gmsh only
    bool clcurv = false;
    std::optional<double> min_size;
    std::optional<double> max_size;
    std::string meshing_algorithm;            // e.g. "del2d", "front2d", "delquad", "del3d", "front3d"
    std::map<std::string, double> gmsh_options; // additional Gmsh options, e.g. {"Mesh.Smoothing", 5}

    Mesh create() const;

    // The backend create() uses
    MeshBackend resolvedBackend() const;

private:
    const Geometry &m_geometry;
};

// Python-style shortcut:
//
//   coords, edof, dofs, bdofs, elementmarkers = create_mesh(g, el_type, el_size_factor, dofs_per_node)
//   auto [coords, edof, dofs, bdofs, element_markers] = create_mesh(g, 2, 0.1, 2);
std::tuple<Matrix, IntMatrix, IntMatrix, BoundaryDofs, std::vector<int>> create_mesh(const Geometry &geometry,
                                                                                      int el_type = 2,
                                                                                      double el_size_factor = 1.0,
                                                                                      int dofs_per_node = 1);

// The .geo description of the geometry as written for Gmsh (also useful for inspecting a
// geometry in the Gmsh GUI). el_type decides whether quadrilaterals are requested.
std::string geo_string(const Geometry &geometry, int el_type = 2);

// Traction on the linear boundary (line) elements of a marker, q = [qx, qy] per unit length
// (calfem-python utils.apply_traction_linear_element). Requires the boundary elements of the mesh.
void apply_traction_linear_element(const BoundaryElements &boundaryElements, const Matrix &coords,
                                   const IntMatrix &dofs, ColVec &F, int marker, const RowVec &q);

} // namespace calfem
