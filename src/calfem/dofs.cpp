//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/dofs.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace calfem {

std::vector<Eigen::Index> toZeroIndex(const std::vector<Eigen::Index> &vec)
{
    std::vector<Eigen::Index> idx;
    idx.reserve(vec.size());

    for (auto &dof : vec)
        idx.push_back(dof - 1);

    return idx;
}

void extractAllDofs(Eigen::Index nDofs, const std::vector<Eigen::Index> &bcDofs, std::vector<Eigen::Index> &freeDofs,
                    std::vector<Eigen::Index> &allDofs)
{
    allDofs.clear();
    freeDofs.clear();

    for (auto i = 1; i <= nDofs; i++)
        allDofs.push_back(i);

    // set_difference requires sorted ranges

    std::vector<Eigen::Index> sortedBcDofs(bcDofs);
    std::sort(sortedBcDofs.begin(), sortedBcDofs.end());
    sortedBcDofs.erase(std::unique(sortedBcDofs.begin(), sortedBcDofs.end()), sortedBcDofs.end());

    std::set_difference(allDofs.begin(), allDofs.end(), sortedBcDofs.begin(), sortedBcDofs.end(),
                        std::back_inserter(freeDofs));
}

void extractElDisp(const IntMatrix &edof, const ColVec &a, Matrix &ed)
{
    auto nElements = edof.rows();
    auto nDofs = edof.cols();

    ed.resize(nElements, nDofs);

    for (auto i = 0; i < nElements; i++)
    {
        IntColVec ind = edof.row(i).array() - 1;
        ed.row(i) = a(ind);
    }
}

void extractEdof(const IntMatrix &dofs, const Matrix &coords, const IntMatrix &topo, IntMatrix &edof, Matrix &ex,
                 Matrix &ey, Matrix &ez)
{
    auto dofsPerNode = dofs.cols();
    auto nElements = topo.rows();
    auto nodesPerElement = topo.cols();
    auto nDims = coords.cols();

    edof.resize(nElements, nodesPerElement * dofsPerNode);
    edof.setZero();

    ex.resize(nElements, nodesPerElement);
    ex.setZero();

    ey.resize(nElements, nodesPerElement);
    ey.setZero();

    ez.resize(nElements, nodesPerElement);
    ez.setZero();

    for (auto i = 0; i < nElements; i++)
    {
        for (auto j = 0; j < nodesPerElement; j++)
        {
            auto k = topo(i, j);

            ex(i, j) = coords(k - 1, 0);

            if (nDims > 1)
                ey(i, j) = coords(k - 1, 1);

            if (nDims > 2)
                ez(i, j) = coords(k - 1, 2);

            edof.block(i, j * dofsPerNode, 1, dofsPerNode) = dofs.row(k - 1);
        }
    }
}

// ---- Value-returning versions ----

std::tuple<std::vector<Eigen::Index>, std::vector<Eigen::Index>> extractAllDofs(
    Eigen::Index nDofs, const std::vector<Eigen::Index> &bcDofs)
{
    std::vector<Eigen::Index> freeDofs, allDofs;
    extractAllDofs(nDofs, bcDofs, freeDofs, allDofs);
    return {freeDofs, allDofs};
}

Matrix extractElDisp(const IntMatrix &edof, const ColVec &a)
{
    Matrix ed;
    extractElDisp(edof, a, ed);
    return ed;
}

std::tuple<IntMatrix, Matrix, Matrix, Matrix> extractEdof(const IntMatrix &dofs, const Matrix &coords,
                                                          const IntMatrix &topo)
{
    IntMatrix edof;
    Matrix ex, ey, ez;
    extractEdof(dofs, coords, topo, edof, ex, ey, ez);
    return {edof, ex, ey, ez};
}

Matrix stress2nodal(const ColVec &eseff, const IntMatrix &edof)
{
    ColVec values = ColVec::Zero(edof.maxCoeff());
    auto elnodes = edof.cols() / 2;

    for (Eigen::Index i = 0; i < edof.rows(); i++)
        for (Eigen::Index j = 0; j < edof.cols(); j++)
            values(edof(i, j) - 1) += eseff(i) / static_cast<double>(elnodes);

    Matrix ev(edof.rows(), elnodes);
    for (Eigen::Index i = 0; i < edof.rows(); i++)
        for (Eigen::Index k = 0; k < elnodes; k++)
            ev(i, k) = values(edof(i, 2 * k) - 1);

    return ev;
}

// ---- Dof numbering and element coordinates ----

IntMatrix create_dofs(int nCoords, int nDof)
{
    IntMatrix dofs(nCoords, nDof);
    for (int i = 0; i < nCoords; i++)
        for (int j = 0; j < nDof; j++)
            dofs(i, j) = i * nDof + j + 1;
    return dofs;
}

std::tuple<Matrix, Matrix, Matrix> coord_extract(const IntMatrix &edof, const Matrix &coords, const IntMatrix &dofs,
                                                 int nen)
{
    auto nElements = edof.rows();
    auto nElementNodes = nen < 0 ? edof.cols() / dofs.cols() : static_cast<Eigen::Index>(nen);
    auto dofsPerNode = edof.cols() / nElementNodes; // dofs per node in edof

    if ((nElementNodes * dofs.cols() != edof.cols()) || (dofsPerNode > dofs.cols()))
        throw std::invalid_argument("coord_extract: edof does not match dofs and the number of element nodes.");

    // Nodes are identified by their first dofsPerNode dofs
    std::map<std::vector<int>, Eigen::Index> nodeIndex;
    for (Eigen::Index i = 0; i < dofs.rows(); i++)
    {
        std::vector<int> key;
        for (Eigen::Index k = 0; k < dofsPerNode; k++)
            key.push_back(dofs(i, k));
        nodeIndex.emplace(key, i);
    }

    Matrix ex = Matrix::Zero(nElements, nElementNodes);
    Matrix ey = Matrix::Zero(nElements, nElementNodes);
    Matrix ez = Matrix::Zero(nElements, nElementNodes);

    for (Eigen::Index e = 0; e < nElements; e++)
        for (Eigen::Index n = 0; n < nElementNodes; n++)
        {
            std::vector<int> key;
            for (Eigen::Index k = 0; k < dofsPerNode; k++)
                key.push_back(edof(e, n * dofsPerNode + k));

            auto it = nodeIndex.find(key);
            if (it == nodeIndex.end())
                throw std::invalid_argument("coord_extract: element " + std::to_string(e + 1) +
                                            " has a node that is not in dofs.");

            ex(e, n) = coords(it->second, 0);
            if (coords.cols() > 1)
                ey(e, n) = coords(it->second, 1);
            if (coords.cols() > 2)
                ez(e, n) = coords(it->second, 2);
        }

    return {ex, ey, ez};
}

std::tuple<Matrix, IntMatrix, IntMatrix> convert_to_node_topo(const IntMatrix &edof, const Matrix &ex,
                                                              const Matrix &ey, const Matrix &ez, int n_dofs_per_node,
                                                              bool ignore_first)
{
    Eigen::Index first = ignore_first ? 1 : 0;
    auto nNodes = (edof.cols() - first) / n_dofs_per_node;

    std::map<std::vector<int>, Eigen::Index> nodeNumber;
    std::vector<std::vector<int>> nodeDofList;
    std::vector<RowVec3> nodeCoords;

    IntMatrix topo(edof.rows(), nNodes);

    for (Eigen::Index e = 0; e < edof.rows(); e++)
        for (Eigen::Index n = 0; n < nNodes; n++)
        {
            std::vector<int> key;
            for (Eigen::Index k = 0; k < n_dofs_per_node; k++)
                key.push_back(edof(e, first + n * n_dofs_per_node + k));

            RowVec3 c(ex(e, n), ey(e, n), ez(e, n));

            auto [it, inserted] = nodeNumber.emplace(key, static_cast<Eigen::Index>(nodeDofList.size()));
            if (inserted)
            {
                nodeDofList.push_back(key);
                nodeCoords.push_back(c);
            }
            else
                nodeCoords[it->second] = c; // last occurrence, as calfem-python

            topo(e, n) = static_cast<int>(it->second);
        }

    auto count = static_cast<Eigen::Index>(nodeDofList.size());
    Matrix coords(count, 3);
    IntMatrix nodeDofs(count, n_dofs_per_node);

    for (Eigen::Index i = 0; i < count; i++)
    {
        coords.row(i) = nodeCoords[i];
        for (Eigen::Index k = 0; k < n_dofs_per_node; k++)
            nodeDofs(i, k) = nodeDofList[i][k];
    }

    return {coords, topo, nodeDofs};
}

} // namespace calfem
