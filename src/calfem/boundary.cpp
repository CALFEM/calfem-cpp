//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Boundary conditions and loads, ported from calfem-python (utils.py).

#include <calfem/boundary.h>

#include <stdexcept>
#include <string>

namespace calfem {

namespace {

// Dofs of the marker; dimension 0 selects all, otherwise every nDim-th dof starting at
// dimension - 1 (as calfem-python's dofs[dimension-1::nDim]).
std::vector<int> selectDofs(const BoundaryDofs &boundaryDofs, int marker, int dimension, int nDim, const char *caller)
{
    auto it = boundaryDofs.find(marker);
    if (it == boundaryDofs.end())
        throw std::invalid_argument(std::string(caller) + ": boundary marker " + std::to_string(marker) +
                                    " does not exist.");

    if ((dimension < 0) || (dimension > nDim))
        throw std::invalid_argument(std::string(caller) + ": invalid dimension " + std::to_string(dimension) + ".");

    if (dimension == 0)
        return it->second;

    std::vector<int> selected;
    for (size_t i = static_cast<size_t>(dimension - 1); i < it->second.size(); i += static_cast<size_t>(nDim))
        selected.push_back(it->second[i]);
    return selected;
}

std::tuple<IntColVec, ColVec> addPrescribed(const IntColVec &bcPrescr, const ColVec &bcVal,
                                            const std::vector<int> &add, double value)
{
    // Sorted, unique; the first occurrence (existing values first) wins, as np.unique
    std::map<int, double> bc;
    for (Eigen::Index i = 0; i < bcPrescr.size(); i++)
        bc.emplace(bcPrescr(i), bcVal(i));
    for (int dof : add)
        bc.emplace(dof, value);

    IntColVec dofs(static_cast<Eigen::Index>(bc.size()));
    ColVec vals(static_cast<Eigen::Index>(bc.size()));
    Eigen::Index i = 0;
    for (const auto &[dof, val] : bc)
    {
        dofs(i) = dof;
        vals(i++) = val;
    }

    return {dofs, vals};
}

void addForce(ColVec &f, const std::vector<int> &dofs, double value)
{
    for (int dof : dofs)
        f(dof - 1) += value;
}

std::vector<int> nodeDofs(int nodeIdx, const IntMatrix &dofs, int dimension, const char *caller)
{
    if ((nodeIdx < 0) || (nodeIdx >= dofs.rows()))
        throw std::invalid_argument(std::string(caller) + ": node index " + std::to_string(nodeIdx) +
                                    " is outside the dofs matrix.");

    if ((dimension < 0) || (dimension > dofs.cols()))
        throw std::invalid_argument(std::string(caller) + ": invalid dimension " + std::to_string(dimension) + ".");

    std::vector<int> selected;
    if (dimension == 0)
        for (Eigen::Index j = 0; j < dofs.cols(); j++)
            selected.push_back(dofs(nodeIdx, j));
    else
        selected.push_back(dofs(nodeIdx, dimension - 1));
    return selected;
}

} // namespace

std::tuple<IntColVec, ColVec> apply_bc(const BoundaryDofs &boundaryDofs, const IntColVec &bcPrescr,
                                       const ColVec &bcVal, int marker, double value, int dimension)
{
    return addPrescribed(bcPrescr, bcVal, selectDofs(boundaryDofs, marker, dimension, 2, "apply_bc"), value);
}

std::tuple<IntColVec, ColVec> apply_bc_3d(const BoundaryDofs &boundaryDofs, const IntColVec &bcPrescr,
                                          const ColVec &bcVal, int marker, double value, int dimension)
{
    return addPrescribed(bcPrescr, bcVal, selectDofs(boundaryDofs, marker, dimension, 3, "apply_bc_3d"), value);
}

std::tuple<IntColVec, ColVec> apply_bc_node(int nodeIdx, const IntMatrix &dofs, const IntColVec &bcPrescr,
                                            const ColVec &bcVal, double value, int dimension)
{
    auto add = nodeDofs(nodeIdx, dofs, dimension, "apply_bc_node");
    auto n = static_cast<Eigen::Index>(add.size());

    IntColVec newDofs(bcPrescr.size() + n);
    ColVec newVals(bcVal.size() + n);
    newDofs << bcPrescr, Eigen::Map<IntColVec>(add.data(), n);
    newVals << bcVal, ColVec::Constant(n, value);

    return {newDofs, newVals};
}

void apply_force(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value, int dimension)
{
    addForce(f, selectDofs(boundaryDofs, marker, dimension, 2, "apply_force"), value);
}

void apply_force_3d(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value, int dimension)
{
    addForce(f, selectDofs(boundaryDofs, marker, dimension, 3, "apply_force_3d"), value);
}

void apply_force_total(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value, int dimension)
{
    auto dofs = selectDofs(boundaryDofs, marker, dimension, 2, "apply_force_total");
    addForce(f, dofs, value / static_cast<double>(dofs.size()));
}

void apply_force_total_3d(const BoundaryDofs &boundaryDofs, ColVec &f, int marker, double value, int dimension)
{
    auto dofs = selectDofs(boundaryDofs, marker, dimension, 3, "apply_force_total_3d");
    addForce(f, dofs, value / static_cast<double>(dofs.size()));
}

void apply_force_node(int nodeIdx, const IntMatrix &dofs, ColVec &f, double value, int dimension)
{
    addForce(f, nodeDofs(nodeIdx, dofs, dimension, "apply_force_node"), value);
}

} // namespace calfem
