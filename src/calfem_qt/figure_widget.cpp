//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include "figure_widget.h"
#include "painter.h"

#include <QAction>
#include <QFileDialog>
#include <QKeyEvent>
#include <QLabel>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPainter>
#include <QStatusBar>
#include <QStyle>
#include <QToolBar>
#include <QWheelEvent>

#include <cmath>
#include <stdexcept>

namespace calfem::vis::qt {

// ---- FigureWidget ----

FigureWidget::FigureWidget(Figure &fig, QWidget *parent) : QWidget(parent), m_figure(fig)
{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    setAttribute(Qt::WA_OpaquePaintEvent);
    setMinimumSize(200, 150);
}

Figure &FigureWidget::figure()
{
    return m_figure;
}

Layout FigureWidget::currentLayout() const
{
    if (m_view)
        return layout(m_figure, width(), height(), *m_view);
    return layout(m_figure, width(), height());
}

void FigureWidget::fit()
{
    m_view.reset();
    update();
}

void FigureWidget::zoom(double factor, const QPointF &center)
{
    Layout l = currentLayout();
    double cx = l.mx(center.x()), cy = l.my(center.y());
    const Bounds &v = l.view;
    m_view = Bounds{cx + (v.xmin - cx) / factor, cx + (v.xmax - cx) / factor, cy + (v.ymin - cy) / factor,
                    cy + (v.ymax - cy) / factor};
    update();
}

void FigureWidget::zoom(double factor)
{
    Layout l = currentLayout();
    zoom(factor, QPointF(l.plotLeft + 0.5 * l.plotWidth, l.plotTop + 0.5 * l.plotHeight));
}

QSize FigureWidget::sizeHint() const
{
    return QSize(int(m_figure.width()), int(m_figure.height()));
}

void FigureWidget::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    paintFigure(p, m_figure, currentLayout(), width(), height());
}

void FigureWidget::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton || event->button() == Qt::MiddleButton)
    {
        m_panning = true;
        m_panStart = event->position();
        m_panView = currentLayout().view;
        setCursor(Qt::ClosedHandCursor);
    }
}

void FigureWidget::mouseMoveEvent(QMouseEvent *event)
{
    Layout l = currentLayout();

    if (m_panning)
    {
        // Model units per pixel from the view at the start of the drag
        double sx = (m_panView.xmax - m_panView.xmin) / l.plotWidth;
        double sy = (m_panView.ymax - m_panView.ymin) / l.plotHeight;
        QPointF d = event->position() - m_panStart;
        m_view = Bounds{m_panView.xmin - d.x() * sx, m_panView.xmax - d.x() * sx, m_panView.ymin + d.y() * sy,
                        m_panView.ymax + d.y() * sy};
        update();
        l = currentLayout();
    }

    if (onCursorMoved)
        onCursorMoved(l.mx(event->position().x()), l.my(event->position().y()));
}

void FigureWidget::mouseReleaseEvent(QMouseEvent *)
{
    m_panning = false;
    unsetCursor();
}

void FigureWidget::mouseDoubleClickEvent(QMouseEvent *)
{
    fit();
}

void FigureWidget::wheelEvent(QWheelEvent *event)
{
    double steps = event->angleDelta().y() / 120.0;
    if (steps != 0.0)
        zoom(std::pow(1.2, steps), event->position());
    event->accept();
}

void FigureWidget::keyPressEvent(QKeyEvent *event)
{
    switch (event->key())
    {
    case Qt::Key_Home:
        fit();
        break;
    case Qt::Key_Plus:
    case Qt::Key_Equal:
        zoom(1.25);
        break;
    case Qt::Key_Minus:
        zoom(1.0 / 1.25);
        break;
    case Qt::Key_G:
        m_figure.set_grid(!m_figure.grid());
        update();
        break;
    default:
        QWidget::keyPressEvent(event);
    }
}

// ---- FigureWindow ----

FigureWindow::FigureWindow(Figure &fig, QWidget *parent) : QMainWindow(parent)
{
    QString title = QString("Figure %1").arg(fig.number());
    if (!fig.title().empty())
        title += ": " + QString::fromStdString(fig.title());
    setWindowTitle(title);

    m_widget = new FigureWidget(fig, this);
    setCentralWidget(m_widget);

    auto toolbar = addToolBar("Figure");
    toolbar->setMovable(false);

    auto fitAction = toolbar->addAction(style()->standardIcon(QStyle::SP_BrowserReload), "Fit");
    fitAction->setToolTip("Fit the view to the figure (Home, double-click)");
    connect(fitAction, &QAction::triggered, this, [this] { m_widget->fit(); });

    auto zoomIn = toolbar->addAction("Zoom in");
    zoomIn->setToolTip("Zoom in (+, mouse wheel)");
    connect(zoomIn, &QAction::triggered, this, [this] { m_widget->zoom(1.25); });

    auto zoomOut = toolbar->addAction("Zoom out");
    zoomOut->setToolTip("Zoom out (-, mouse wheel)");
    connect(zoomOut, &QAction::triggered, this, [this] { m_widget->zoom(1.0 / 1.25); });

    auto grid = toolbar->addAction("Grid");
    grid->setCheckable(true);
    grid->setChecked(fig.grid());
    grid->setToolTip("Show the grid (G)");
    connect(grid, &QAction::toggled, this, [this](bool on) {
        m_widget->figure().set_grid(on);
        m_widget->update();
    });

    toolbar->addSeparator();
    auto save = toolbar->addAction(style()->standardIcon(QStyle::SP_DialogSaveButton), "Save");
    save->setShortcut(QKeySequence::Save);
    save->setToolTip("Save the view as SVG, PNG or PDF (Ctrl+S)");
    connect(save, &QAction::triggered, this, [this] { saveAs(); });

    m_position = new QLabel(this);
    statusBar()->addWidget(m_position);
    m_widget->onCursorMoved = [this](double x, double y) {
        m_position->setText(QString("x = %1, y = %2").arg(x, 0, 'g', 6).arg(y, 0, 'g', 6));
    };

    m_widget->setFocus();
}

FigureWidget *FigureWindow::figureWidget() const
{
    return m_widget;
}

void FigureWindow::saveAs()
{
    QString filename =
        QFileDialog::getSaveFileName(this, "Save figure", QString("figure%1.svg").arg(m_widget->figure().number()),
                                     "SVG image (*.svg);;PNG image (*.png);;PDF document (*.pdf);;JPEG image (*.jpg)");
    if (filename.isEmpty())
        return;

    try
    {
        // Saves the visible region at the widget size
        Layout l = m_widget->currentLayout();
        if (!saveFigure(m_widget->figure(), filename, m_widget->width(), m_widget->height(), l.view))
            QMessageBox::warning(this, "Save figure", "Unsupported file format: " + filename);
    }
    catch (const std::exception &e)
    {
        QMessageBox::warning(this, "Save figure", QString::fromStdString(e.what()));
    }
}

} // namespace calfem::vis::qt
