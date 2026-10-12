//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Gmsh backend: meshes the geometry with the Gmsh library (in-process, no executable).
//
// NOTE: Gmsh is licensed under the GPL (version 2 or later). Programs linked with calfem built
// with CALFEM_WITH_GMSH must comply with the GPL when distributed.

#include "mesh_internal.h"

#include <stdexcept>
#include <string>

#ifdef CALFEM_WITH_GMSH

#ifdef CALFEM_GMSH_CWRAP
#include <gmsh.h_cwrap> // C++ API on top of the C API (compiler independent)
#else
#include <gmsh.h>
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace calfem::detail {

namespace {

const std::set<int> secondOrderElements = {8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19};

const std::map<std::string, int> algorithms2D = {{"meshadapt", 1}, {"auto", 2},    {"initial2d", 3},
                                                 {"del2d", 5},     {"front2d", 6}, {"delquad", 8},
                                                 {"quadqs", 11}};
const std::map<std::string, int> algorithms3D = {
    {"del3d", 1}, {"initial3d", 3}, {"front3d", 4}, {"mmg3d", 7}, {"hxt", 10}};

// Initializes Gmsh unless the application already has, and finalizes it again
class GmshSession {
public:
    GmshSession() : m_owner(!gmsh::isInitialized())
    {
        if (m_owner)
        {
            // Gmsh expects at least the program name (as calfem-python's sys.argv[:1])
            static char name[] = "calfem";
            static char *argv[] = {name, nullptr};
            gmsh::initialize(1, argv, false, false);
        }
    }

    ~GmshSession()
    {
        if (m_owner)
            gmsh::finalize();
    }

    GmshSession(const GmshSession &) = delete;
    GmshSession &operator=(const GmshSession &) = delete;

private:
    bool m_owner;
};

// Temporary .geo file, removed when done
class TempGeoFile {
public:
    explicit TempGeoFile(const std::string &contents)
    {
        auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        m_path = std::filesystem::temp_directory_path() / ("calfem_" + std::to_string(stamp) + ".geo");

        std::ofstream f(m_path);
        if (!f)
            throw std::runtime_error("Mesh: could not write the temporary geometry file " + m_path.string() + ".");
        f << contents;
    }

    ~TempGeoFile()
    {
        std::error_code ec;
        std::filesystem::remove(m_path, ec);
    }

    std::string path() const
    {
        return m_path.string();
    }

private:
    std::filesystem::path m_path;
};

} // namespace

RawMesh meshGmsh(const Geometry &geometry, const MeshGenerator &gen)
{
    int dim = geometry.is3D() ? 3 : 2;
    TempGeoFile geo(geo_string(geometry, gen.el_type));

    GmshSession session;

    try
    {
        gmsh::option::setNumber("General.Verbosity", 0);
        gmsh::option::setNumber("General.Terminal", 0);

        gmsh::open(geo.path());
        gmsh::model::geo::synchronize();

        if (secondOrderElements.count(gen.el_type))
            gmsh::option::setNumber("Mesh.ElementOrder", 2);

        if (!gen.meshing_algorithm.empty())
        {
            std::string algorithm = gen.meshing_algorithm;
            std::transform(algorithm.begin(), algorithm.end(), algorithm.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

            if (algorithms3D.count(algorithm))
                gmsh::option::setNumber("Mesh.Algorithm3D", algorithms3D.at(algorithm));
            else if (algorithms2D.count(algorithm))
                gmsh::option::setNumber("Mesh.Algorithm", algorithms2D.at(algorithm));
            else
                throw std::invalid_argument("Mesh: unknown meshing_algorithm '" + gen.meshing_algorithm + "'.");
        }

        gmsh::option::setNumber("Mesh.MeshSizeFactor", gen.el_size_factor);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", gen.clcurv ? 1 : 0);
        if (gen.min_size)
            gmsh::option::setNumber("Mesh.MeshSizeMin", *gen.min_size);
        if (gen.max_size)
            gmsh::option::setNumber("Mesh.MeshSizeMax", *gen.max_size);
        for (const auto &[name, value] : gen.gmsh_options)
            gmsh::option::setNumber(name, value);

        gmsh::model::mesh::generate(dim);

        RawMesh raw;
        raw.dim = dim;

        // ---- Elements of the physical groups ----
        //
        // As calfem-python: elements of type el_type are mesh elements (marker = physical group),
        // all other elements contribute their nodes to bdofs of their physical group. Only nodes
        // used by these elements are kept (as in the .msh file calfem-python reads), numbered in
        // order of their Gmsh tags.

        struct GroupElement {
            int type;
            std::size_t tag;
            std::vector<std::size_t> nodeTags;
            int marker;
        };

        std::vector<GroupElement> groupElements;
        std::set<std::size_t> usedNodes;

        gmsh::vectorpair groups;
        gmsh::model::getPhysicalGroups(groups);

        for (const auto &[gdim, marker] : groups)
        {
            std::vector<int> entities;
            gmsh::model::getEntitiesForPhysicalGroup(gdim, marker, entities);

            for (int entity : entities)
            {
                std::vector<int> types;
                std::vector<std::vector<std::size_t>> elementTags, elementNodes;
                gmsh::model::mesh::getElements(types, elementTags, elementNodes, gdim, entity);

                for (size_t t = 0; t < types.size(); t++)
                {
                    auto nElements = elementTags[t].size();
                    if (nElements == 0)
                        continue;
                    auto nn = elementNodes[t].size() / nElements;

                    for (size_t e = 0; e < nElements; e++)
                    {
                        std::vector<std::size_t> nodes(elementNodes[t].begin() + e * nn,
                                                       elementNodes[t].begin() + (e + 1) * nn);
                        usedNodes.insert(nodes.begin(), nodes.end());
                        groupElements.push_back({types[t], elementTags[t][e], nodes, marker});
                    }
                }
            }
        }

        // ---- Nodes ----

        std::vector<std::size_t> nodeTags;
        std::vector<double> coord, parametric;
        gmsh::model::mesh::getNodes(nodeTags, coord, parametric, -1, -1, false, false);

        std::map<std::size_t, size_t> tagIndex;
        for (size_t i = 0; i < nodeTags.size(); i++)
            tagIndex[nodeTags[i]] = i;

        std::map<std::size_t, int> nodeNumber; // gmsh tag -> 1-based node number
        for (std::size_t tag : usedNodes) // ascending tags
        {
            size_t i = tagIndex.at(tag);
            nodeNumber[tag] = static_cast<int>(raw.nodes.size()) + 1;
            raw.nodes.emplace_back(coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);
        }

        std::vector<std::tuple<std::size_t, std::vector<int>, int>> elements; // tag, nodes, marker

        for (const auto &ge : groupElements)
        {
            std::vector<int> nodes;
            for (std::size_t tag : ge.nodeTags)
                nodes.push_back(nodeNumber.at(tag));

            if (ge.type == gen.el_type)
                elements.emplace_back(ge.tag, nodes, ge.marker);
            else
            {
                raw.boundaryNodes[ge.marker].insert(nodes.begin(), nodes.end());
                raw.boundaryElements[ge.marker].push_back({ge.type, nodes});
            }
        }

        std::sort(elements.begin(), elements.end(),
                  [](const auto &a, const auto &b) { return std::get<0>(a) < std::get<0>(b); });

        for (const auto &[tag, nodes, marker] : elements)
        {
            raw.elements.push_back(nodes);
            raw.markers.push_back(marker);
        }

        if (raw.elements.empty())
            throw std::runtime_error("Mesh: Gmsh produced no elements of type " + std::to_string(gen.el_type) + ".");

        return raw;
    }
    catch (const std::invalid_argument &)
    {
        throw;
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(std::string("Mesh: Gmsh failed: ") + e.what());
    }
}

} // namespace calfem::detail

#else

namespace calfem::detail {

RawMesh meshGmsh(const Geometry &, const MeshGenerator &)
{
    throw std::runtime_error("Mesh: calfem was built without Gmsh support (CALFEM_WITH_GMSH). Use the built-in "
                             "mesher (MeshBackend::BuiltIn) or rebuild with the Gmsh SDK.");
}

} // namespace calfem::detail

#endif
