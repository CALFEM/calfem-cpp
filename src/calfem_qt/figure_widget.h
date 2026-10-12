//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Interactive figure widget and window (private to calfem::vis_qt).

#pragma once

#include <calfem/vis.h>

#include <QMainWindow>
#include <QPointF>
#include <QWidget>

#include <functional>
#include <optional>

class QLabel;

namespace calfem::vis::qt {

// Draws a Figure; drag with the left or middle mouse button to pan, wheel to zoom around the
// cursor, double-click or Home to fit the view to the figure, +/- to zoom.
class FigureWidget : public QWidget {
public:
    explicit FigureWidget(Figure &fig, QWidget *parent = nullptr);

    Figure &figure();

    // Current layout (visible model region, ticks) for the widget size.
    Layout currentLayout() const;

    // Back to the default view (axis limits or data bounds).
    void fit();
    // Zooms by factor (> 1 zooms in) around the pixel position center.
    void zoom(double factor, const QPointF &center);
    void zoom(double factor);

    // Called with the model coordinates under the cursor.
    std::function<void(double x, double y)> onCursorMoved;

    QSize sizeHint() const override;

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

private:
    Figure &m_figure;
    std::optional<Bounds> m_view;
    bool m_panning = false;
    QPointF m_panStart;
    Bounds m_panView;
};

// Main window with a FigureWidget, a toolbar (fit, zoom, grid, save) and the cursor position in
// the status bar.
class FigureWindow : public QMainWindow {
public:
    explicit FigureWindow(Figure &fig, QWidget *parent = nullptr);

    FigureWidget *figureWidget() const;

private:
    void saveAs();

    FigureWidget *m_widget;
    QLabel *m_position;
};

} // namespace calfem::vis::qt
