//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Drawing of calfem::vis figures with QPainter (private to calfem::vis_qt).

#pragma once

#include <calfem/vis.h>

#include <QString>

class QPainter;

namespace calfem::vis::qt {

// Draws fig in the painter's width x height pixel area using layout l (axes, ticks, grid,
// items and title).
void paintFigure(QPainter &painter, const Figure &fig, const Layout &l, double width, double height);

// Writes fig (width x height pixels, model region view) to filename. Supports .svg, .png,
// .jpg/.jpeg, .bmp and .pdf; returns false for other extensions. Throws std::runtime_error when
// the file cannot be written.
bool saveFigure(const Figure &fig, const QString &filename, double width, double height, const Bounds &view);
bool saveFigure(const Figure &fig, const QString &filename);

// Creates the QApplication if the program has none (needed for windows, fonts and images).
void ensureApplication();

} // namespace calfem::vis::qt
