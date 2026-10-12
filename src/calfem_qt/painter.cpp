//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include "painter.h"

#include <QApplication>
#include <QFileInfo>
#include <QFont>
#include <QFontMetricsF>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QPdfWriter>
#include <QPolygonF>

#include <cmath>
#include <fstream>
#include <stdexcept>

namespace calfem::vis::qt {

namespace {

const double pi = 3.14159265358979323846;

QColor toQColor(const Color &c)
{
    return QColor::fromRgbF(float(c.r), float(c.g), float(c.b));
}

QPen linePen(const Line &line)
{
    QPen pen(toQColor(line.color), line.width);
    pen.setJoinStyle(Qt::RoundJoin);
    pen.setCapStyle(Qt::FlatCap);
    // Dash patterns are in units of the pen width (as the SVG output)
    if (line.style == LineStyle::Dashed)
        pen.setDashPattern({5.0, 5.0});
    else if (line.style == LineStyle::Dotted)
        pen.setDashPattern({1.0, 1.65});
    return pen;
}

QFont font(double pointSize)
{
    QFont f;
    f.setPixelSize(std::max(1, int(std::lround(pointSize))));
    return f;
}

void paintLine(QPainter &p, const Line &line, const Layout &l)
{
    if (line.points.size() < 2)
        return;
    QPolygonF poly;
    for (const auto &pt : line.points)
        poly << QPointF(l.px(pt.x), l.py(pt.y));

    p.setPen(linePen(line));
    p.setBrush(Qt::NoBrush);
    if (line.closed)
        p.drawPolygon(poly);
    else
        p.drawPolyline(poly);
}

void paintMarkers(QPainter &p, const Markers &markers, const Layout &l)
{
    if (markers.marker == Marker::None)
        return;

    QColor color = toQColor(markers.color);
    p.setPen(QPen(color, 1.0));
    p.setBrush(markers.filled ? QBrush(color) : Qt::NoBrush);
    double r = 0.5 * markers.size;

    for (const auto &pt : markers.points)
    {
        QPointF c(l.px(pt.x), l.py(pt.y));
        switch (markers.marker)
        {
        case Marker::Circle:
            p.drawEllipse(c, r, r);
            break;
        case Marker::Dot:
            p.setBrush(color);
            p.drawEllipse(c, 0.3 * r, 0.3 * r);
            p.setBrush(markers.filled ? QBrush(color) : Qt::NoBrush);
            break;
        case Marker::Cross:
            p.drawLine(c + QPointF(-r, -r), c + QPointF(r, r));
            p.drawLine(c + QPointF(-r, r), c + QPointF(r, -r));
            break;
        case Marker::Star: {
            QPolygonF star;
            for (int k = 0; k < 10; k++)
            {
                double angle = -0.5 * pi + k * pi / 5.0;
                double radius = k % 2 == 0 ? r : 0.4 * r;
                star << c + QPointF(radius * std::cos(angle), radius * std::sin(angle));
            }
            p.drawPolygon(star);
            break;
        }
        default:
            break;
        }
    }
}

// Draws text with its anchor at (x, y) in pixels
void paintText(QPainter &p, const QString &text, double x, double y, double fontSize, HAlign halign, VAlign valign)
{
    QFont f = font(fontSize);
    QFontMetricsF fm(f);
    double w = fm.horizontalAdvance(text);

    if (halign == HAlign::Center)
        x -= 0.5 * w;
    else if (halign == HAlign::Right)
        x -= w;

    switch (valign)
    {
    case VAlign::Bottom:
        y -= fm.descent();
        break;
    case VAlign::Center:
        y += 0.5 * (fm.ascent() - fm.descent());
        break;
    case VAlign::Top:
        y += fm.ascent();
        break;
    default:
        break;
    }

    p.setFont(f);
    p.drawText(QPointF(x, y), text);
}

void paintAxes(QPainter &p, const Layout &l, bool grid)
{
    QRectF plot(l.plotLeft, l.plotTop, l.plotWidth, l.plotHeight);
    double bottom = plot.bottom();

    if (grid)
    {
        p.setPen(QPen(QColor(208, 208, 208), 0.8));
        for (double t : l.xticks)
            p.drawLine(QPointF(l.px(t), plot.top()), QPointF(l.px(t), bottom));
        for (double t : l.yticks)
            p.drawLine(QPointF(plot.left(), l.py(t)), QPointF(plot.right(), l.py(t)));
    }

    p.setPen(QPen(Qt::black, 0.8));
    p.setBrush(Qt::NoBrush);
    p.drawRect(plot);

    for (double t : l.xticks)
    {
        p.drawLine(QPointF(l.px(t), bottom), QPointF(l.px(t), bottom + 4.0));
        paintText(p, QString::fromStdString(format_tick(t, l.tickStepX)), l.px(t), bottom + 16.0, 10.0,
                  HAlign::Center, VAlign::Baseline);
    }
    for (double t : l.yticks)
    {
        p.drawLine(QPointF(plot.left() - 4.0, l.py(t)), QPointF(plot.left(), l.py(t)));
        paintText(p, QString::fromStdString(format_tick(t, l.tickStepY)), plot.left() - 7.0, l.py(t), 10.0,
                  HAlign::Right, VAlign::Center);
    }
}

} // namespace

void paintFigure(QPainter &p, const Figure &fig, const Layout &l, double width, double height)
{
    p.save();
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setRenderHint(QPainter::TextAntialiasing, true);
    p.fillRect(QRectF(0.0, 0.0, width, height), Qt::white);

    paintAxes(p, l, fig.grid());

    p.save();
    p.setClipRect(QRectF(l.plotLeft, l.plotTop, l.plotWidth, l.plotHeight));
    for (const auto &item : fig.items())
    {
        if (auto line = std::get_if<Line>(&item))
            paintLine(p, *line, l);
        else if (auto markers = std::get_if<Markers>(&item))
            paintMarkers(p, *markers, l);
        else if (auto text = std::get_if<Text>(&item))
        {
            p.setPen(toQColor(text->color));
            paintText(p, QString::fromStdString(text->text), l.px(text->position.x), l.py(text->position.y),
                      text->fontSize, text->halign, text->valign);
        }
    }
    p.restore();

    if (!fig.title().empty())
    {
        p.setPen(Qt::black);
        paintText(p, QString::fromStdString(fig.title()), l.plotLeft + 0.5 * l.plotWidth, 0.5 * l.plotTop + 6.0,
                  14.0, HAlign::Center, VAlign::Baseline);
    }
    p.restore();
}

void ensureApplication()
{
    if (QCoreApplication::instance() != nullptr)
        return;

    // Kept for the lifetime of the program (destroying a QApplication during static
    // destruction is not safe).
    static int argc = 1;
    static char name[] = "calfem";
    static char *argv[] = {name, nullptr};
    new QApplication(argc, argv);
}

bool saveFigure(const Figure &fig, const QString &filename, double width, double height, const Bounds &view)
{
    QString ext = QFileInfo(filename).suffix().toLower();
    Layout l = layout(fig, width, height, view);

    if (ext == "svg")
    {
        std::ofstream out(filename.toStdString(), std::ios::binary);
        out << to_svg(fig, width, height, view);
        if (!out)
            throw std::runtime_error("save_figure: cannot write " + filename.toStdString() + ".");
        return true;
    }

    if (ext == "png" || ext == "jpg" || ext == "jpeg" || ext == "bmp")
    {
        ensureApplication();
        // Rendered at twice the size for a sharper image
        const double scale = 2.0;
        QImage image(int(std::lround(width * scale)), int(std::lround(height * scale)), QImage::Format_RGB32);
        image.setDevicePixelRatio(scale);
        {
            QPainter p(&image);
            paintFigure(p, fig, l, width, height);
        }
        if (!image.save(filename))
            throw std::runtime_error("save_figure: cannot write " + filename.toStdString() + ".");
        return true;
    }

    if (ext == "pdf")
    {
        ensureApplication();
        QPdfWriter pdf(filename);
        pdf.setResolution(72);
        pdf.setPageSize(QPageSize(QSizeF(width, height), QPageSize::Point));
        pdf.setPageMargins(QMarginsF(0, 0, 0, 0));
        QPainter p;
        if (!p.begin(&pdf))
            throw std::runtime_error("save_figure: cannot write " + filename.toStdString() + ".");
        paintFigure(p, fig, l, width, height);
        p.end();
        return true;
    }

    return false;
}

bool saveFigure(const Figure &fig, const QString &filename)
{
    Layout l = layout(fig, fig.width(), fig.height());
    return saveFigure(fig, filename, fig.width(), fig.height(), l.view);
}

} // namespace calfem::vis::qt
