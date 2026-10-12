//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Qt backend for calfem::vis: show() opens a FigureWindow per figure, save_figure() writes
// raster images and PDF.

#include <calfem/vis_qt.h>

#include "figure_widget.h"
#include "painter.h"

#include <QApplication>
#include <QEventLoop>
#include <QPointer>

#include <memory>
#include <vector>

namespace calfem::vis {

namespace qt {

namespace {

void showFigures(const std::vector<Figure *> &figs, bool block)
{
    ensureApplication();

    std::vector<QPointer<FigureWindow>> windows;
    int offset = 0;
    for (auto fig : figs)
    {
        auto window = new FigureWindow(*fig);
        window->setAttribute(Qt::WA_DeleteOnClose);
        window->resize(window->sizeHint());
        window->move(window->pos() + QPoint(offset, offset));
        window->show();
        windows.push_back(window);
        offset += 30;
    }

    if (!block)
        return;

    // Runs until all windows are closed (also inside a running application event loop)
    QEventLoop loop;
    auto openWindows = std::make_shared<int>(int(windows.size()));
    for (auto &window : windows)
        QObject::connect(window, &QObject::destroyed, &loop, [&loop, openWindows] {
            if (--*openWindows == 0)
                loop.quit();
        });
    loop.exec();
}

bool saveWithQt(const Figure &fig, const std::string &filename)
{
    return saveFigure(fig, QString::fromStdString(filename));
}

} // namespace

void install()
{
    set_backend(Backend{showFigures, saveWithQt});
}

} // namespace qt

QWidget *figure_widget(Figure &fig, QWidget *parent)
{
    qt::ensureApplication();
    return new qt::FigureWidget(fig, parent);
}

} // namespace calfem::vis
