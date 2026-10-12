//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Qt backend for calfem::vis (target calfem::vis_qt, built with CALFEM_WITH_QT).
//
// Including this header installs the Qt backend: vis::show() then opens one window per figure
// (pan with the mouse, zoom with the wheel, double-click or Home to fit, save as SVG, PNG or
// PDF) and save_figure() also writes .png, .jpg, .bmp and .pdf.
//
// show() creates a QApplication if the program has none and must be called from the main
// thread. In a program that runs its own QApplication, figure_widget() embeds a figure in a
// window of the program.

#pragma once

#include <calfem/vis.h>

class QWidget;

namespace calfem::vis {

// Widget drawing fig, as calfem-python's figure_widget(fig, parent). fig must outlive the
// widget; call update() on the widget after changing fig.
QWidget *figure_widget(Figure &fig, QWidget *parent = nullptr);

namespace qt {

// Makes the Qt backend the current backend (done when this header is included).
void install();

namespace detail {
inline const bool installed = (install(), true);
} // namespace detail

} // namespace qt

} // namespace calfem::vis
