//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Matrix and vector types used throughout the library.

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstddef>
#include <vector>

namespace calfem {

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVec;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ColVec;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::MatrixXd Mat;
typedef Eigen::SparseMatrix<double> SpMatrix;

typedef Eigen::Matrix<int, 1, Eigen::Dynamic> IntRowVec;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IntColVec;
typedef Eigen::MatrixXi IntMatrix;
typedef Eigen::MatrixXi IntMat;
typedef Eigen::SparseMatrix<int> SpIntMatrix;

typedef Eigen::Triplet<double, std::ptrdiff_t> Triplet;
typedef std::vector<Triplet> TripletList;

// Fixed-size types for stack allocation (eliminates heap overhead)
typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 12, 12> Matrix12x12;
typedef Eigen::Matrix<double, 12, 1> Vector12;
typedef Eigen::Matrix<double, 6, 12> Matrix6x12;
typedef Eigen::Matrix<double, 6, 1> Vector6;
typedef Eigen::Matrix<double, 4, 12> Matrix4x12;
typedef Eigen::Matrix<double, 4, 1> Vector4;
typedef Eigen::Matrix<double, 1, 3> RowVec3;

} // namespace calfem
