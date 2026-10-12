//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Mesh generation front end, ported from calfem-python (mesh.py, utils.py).

#include <calfem/dofs.h>
#include <calfem/mesh.h>

#include "../detail/log.h"
#include "mesh_internal.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace calfem {

bool gmsh_available()
{
#ifdef CALFEM_WITH_GMSH
    return true;
#else
    return false;
#endif
}

MeshGenerator::MeshGenerator(const Geometry &geometry) : m_geometry(geometry)
{}

MeshBackend MeshGenerator::resolvedBackend() const
{
    if (backend == MeshBackend::Auto)
        return gmsh_available() ? MeshBackend::Gmsh : MeshBackend::BuiltIn;
    return backend;
}

Mesh MeshGenerator::create() const
{
    if (dofs_per_node < 1)
        throw std::invalid_argument("Mesh: dofs_per_node must be at least 1.");
    if (!(el_size_factor > 0.0))
        throw std::invalid_argument("Mesh: el_size_factor must be positive.");

    detail::RawMesh raw;

    if (resolvedBackend() == MeshBackend::Gmsh)
        raw = detail::meshGmsh(m_geometry, *this);
    else
    {
        if (clcurv || min_size || max_size || !meshing_algorithm.empty() || !gmsh_options.empty())
            detail::logMessage("Mesh: warning, clcurv, min_size, max_size, meshing_algorithm and gmsh_options are "
                               "only used by the Gmsh backend.");
        raw = detail::meshBuiltIn(m_geometry, el_type, el_size_factor);
    }

    return detail::finalizeMesh(raw, dofs_per_node);
}

std::tuple<Matrix, IntMatrix, IntMatrix, BoundaryDofs, std::vector<int>> create_mesh(const Geometry &geometry,
                                                                                      int el_type,
                                                                                      double el_size_factor,
                                                                                      int dofs_per_node)
{
    MeshGenerator gen(geometry);
    gen.el_type = el_type;
    gen.el_size_factor = el_size_factor;
    gen.dofs_per_node = dofs_per_node;

    Mesh m = gen.create();
    return {m.coords, m.edof, m.dofs, m.bdofs, m.element_markers};
}

namespace detail {

Mesh finalizeMesh(const RawMesh &raw, int dofsPerNode)
{
    Mesh m;

    auto nNodes = static_cast<Eigen::Index>(raw.nodes.size());
    m.coords.resize(nNodes, raw.dim);
    for (Eigen::Index i = 0; i < nNodes; i++)
        m.coords.row(i) = raw.nodes[i].head(raw.dim).transpose();

    m.dofs = create_dofs(static_cast<int>(nNodes), dofsPerNode);

    auto nElements = static_cast<Eigen::Index>(raw.elements.size());
    auto nodesPerElement = nElements > 0 ? static_cast<Eigen::Index>(raw.elements[0].size()) : 0;

    m.topo.resize(nElements, nodesPerElement);
    m.edof.resize(nElements, nodesPerElement * dofsPerNode);

    for (Eigen::Index e = 0; e < nElements; e++)
        for (Eigen::Index k = 0; k < nodesPerElement; k++)
        {
            int node = raw.elements[e][k];
            m.topo(e, k) = node;
            for (int d = 0; d < dofsPerNode; d++)
                m.edof(e, k * dofsPerNode + d) = m.dofs(node - 1, d);
        }

    m.element_markers = raw.markers;

    for (const auto &[marker, nodes] : raw.boundaryNodes)
        for (int node : nodes)
            for (int d = 0; d < dofsPerNode; d++)
                m.bdofs[marker].push_back(m.dofs(node - 1, d));

    m.boundary_elements = raw.boundaryElements;
    return m;
}

} // namespace detail

// ---- Traction on boundary line elements ----

void apply_traction_linear_element(const BoundaryElements &boundaryElements, const Matrix &coords,
                                   const IntMatrix &dofs, ColVec &F, int marker, const RowVec &q)
{
    auto it = boundaryElements.find(marker);
    if (it == boundaryElements.end())
        throw std::invalid_argument("apply_traction_linear_element: boundary marker " + std::to_string(marker) +
                                    " does not exist.");

    if (dofs.cols() != 2)
        throw std::invalid_argument("apply_traction_linear_element: two dofs per node are required.");

    for (const auto &element : it->second)
        if (element.elm_type != 1)
            throw std::invalid_argument("apply_traction_linear_element: marker " + std::to_string(marker) +
                                        " has elements that are not 2-node lines.");

    const double xi[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    ColVec qv = q.transpose();

    for (const auto &element : it->second)
    {
        int n1 = element.nodes[0], n2 = element.nodes[1];
        double J = (coords.row(n1 - 1) - coords.row(n2 - 1)).norm() / 2.0;

        ColVec f = ColVec::Zero(4);
        for (double x : xi)
        {
            double N1 = 1.0 - (1.0 + x) / 2.0, N2 = (1.0 + x) / 2.0;
            Matrix N(2, 4);
            N << N1, 0.0, N2, 0.0, //
                0.0, N1, 0.0, N2;
            f += N.transpose() * qv * J;
        }

        int idx[4] = {dofs(n1 - 1, 0), dofs(n1 - 1, 1), dofs(n2 - 1, 0), dofs(n2 - 1, 1)};
        for (int k = 0; k < 4; k++)
            F(idx[k] - 1) += f(k);
    }
}

} // namespace calfem
