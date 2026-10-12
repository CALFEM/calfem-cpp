//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// SVG output of figures.

#include <calfem/vis.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <locale>
#include <sstream>

namespace calfem::vis {

namespace {

std::string svgColor(const Color &c)
{
    auto component = [](double v) { return int(std::lround(std::clamp(v, 0.0, 1.0) * 255.0)); };
    char buffer[8];
    std::snprintf(buffer, sizeof(buffer), "#%02x%02x%02x", component(c.r), component(c.g), component(c.b));
    return buffer;
}

std::string escape(const std::string &text)
{
    std::string result;
    for (char c : text)
    {
        switch (c)
        {
        case '&':
            result += "&amp;";
            break;
        case '<':
            result += "&lt;";
            break;
        case '>':
            result += "&gt;";
            break;
        case '"':
            result += "&quot;";
            break;
        default:
            result += c;
        }
    }
    return result;
}

class SvgWriter {
public:
    SvgWriter(const Figure &fig, double width, double height, const Layout &layout)
        : m_fig(fig), m_width(width), m_height(height), m_layout(layout)
    {
        m_out.imbue(std::locale::classic());
        m_out << std::setprecision(6);
    }

    std::string write()
    {
        const auto &l = m_layout;

        m_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        m_out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << m_width << "\" height=\"" << m_height
              << "\" viewBox=\"0 0 " << m_width << " " << m_height << "\" font-family=\"sans-serif\">\n";
        m_out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        m_out << "<defs><clipPath id=\"plot\"><rect x=\"" << l.plotLeft << "\" y=\"" << l.plotTop << "\" width=\""
              << l.plotWidth << "\" height=\"" << l.plotHeight << "\"/></clipPath></defs>\n";

        writeGrid();
        writeAxes();

        m_out << "<g clip-path=\"url(#plot)\" fill=\"none\" stroke-linejoin=\"round\" stroke-linecap=\"butt\">\n";
        for (const auto &item : m_fig.items())
        {
            if (auto line = std::get_if<Line>(&item))
                writeLine(*line);
            else if (auto markers = std::get_if<Markers>(&item))
                writeMarkers(*markers);
            else if (auto text = std::get_if<Text>(&item))
                writeText(*text);
        }
        m_out << "</g>\n";

        if (!m_fig.title().empty())
            m_out << "<text x=\"" << l.plotLeft + 0.5 * l.plotWidth << "\" y=\"" << 0.5 * l.plotTop + 6.0
                  << "\" font-size=\"14\" text-anchor=\"middle\">" << escape(m_fig.title()) << "</text>\n";

        m_out << "</svg>\n";
        return m_out.str();
    }

private:
    void writeGrid()
    {
        if (!m_fig.grid())
            return;
        const auto &l = m_layout;
        m_out << "<g stroke=\"#d0d0d0\" stroke-width=\"0.8\">\n";
        for (double t : l.xticks)
            m_out << "<line x1=\"" << l.px(t) << "\" y1=\"" << l.plotTop << "\" x2=\"" << l.px(t) << "\" y2=\""
                  << l.plotTop + l.plotHeight << "\"/>\n";
        for (double t : l.yticks)
            m_out << "<line x1=\"" << l.plotLeft << "\" y1=\"" << l.py(t) << "\" x2=\"" << l.plotLeft + l.plotWidth
                  << "\" y2=\"" << l.py(t) << "\"/>\n";
        m_out << "</g>\n";
    }

    void writeAxes()
    {
        const auto &l = m_layout;
        double bottom = l.plotTop + l.plotHeight;

        m_out << "<rect x=\"" << l.plotLeft << "\" y=\"" << l.plotTop << "\" width=\"" << l.plotWidth
              << "\" height=\"" << l.plotHeight << "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.8\"/>\n";

        m_out << "<g stroke=\"black\" stroke-width=\"0.8\">\n";
        for (double t : l.xticks)
            m_out << "<line x1=\"" << l.px(t) << "\" y1=\"" << bottom << "\" x2=\"" << l.px(t) << "\" y2=\""
                  << bottom + 4.0 << "\"/>\n";
        for (double t : l.yticks)
            m_out << "<line x1=\"" << l.plotLeft - 4.0 << "\" y1=\"" << l.py(t) << "\" x2=\"" << l.plotLeft
                  << "\" y2=\"" << l.py(t) << "\"/>\n";
        m_out << "</g>\n";

        m_out << "<g font-size=\"10\" fill=\"black\">\n";
        for (double t : l.xticks)
            m_out << "<text x=\"" << l.px(t) << "\" y=\"" << bottom + 16.0 << "\" text-anchor=\"middle\">"
                  << format_tick(t, l.tickStepX) << "</text>\n";
        for (double t : l.yticks)
            m_out << "<text x=\"" << l.plotLeft - 7.0 << "\" y=\"" << l.py(t) + 3.5 << "\" text-anchor=\"end\">"
                  << format_tick(t, l.tickStepY) << "</text>\n";
        m_out << "</g>\n";
    }

    void writePoints(const std::vector<Point> &points)
    {
        for (std::size_t i = 0; i < points.size(); i++)
            m_out << (i > 0 ? " " : "") << m_layout.px(points[i].x) << "," << m_layout.py(points[i].y);
    }

    void writeLine(const Line &line)
    {
        if (line.points.size() < 2)
            return;
        m_out << (line.closed ? "<polygon" : "<polyline") << " points=\"";
        writePoints(line.points);
        m_out << "\" stroke=\"" << svgColor(line.color) << "\" stroke-width=\"" << line.width << "\"";
        writeDashArray(line.style, line.width);
        m_out << "/>\n";
    }

    void writeDashArray(LineStyle style, double width)
    {
        if (style == LineStyle::Dashed)
            m_out << " stroke-dasharray=\"" << 5.0 * width << "," << 5.0 * width << "\"";
        else if (style == LineStyle::Dotted)
            m_out << " stroke-dasharray=\"" << width << "," << 1.65 * width << "\"";
    }

    void writeMarkers(const Markers &markers)
    {
        if (markers.marker == Marker::None || markers.points.empty())
            return;

        std::string color = svgColor(markers.color);
        std::string fill = markers.filled ? color : "none";
        double r = 0.5 * markers.size;

        m_out << "<g stroke=\"" << color << "\" stroke-width=\"1\" fill=\"" << fill << "\">\n";
        for (const auto &p : markers.points)
        {
            double x = m_layout.px(p.x), y = m_layout.py(p.y);
            switch (markers.marker)
            {
            case Marker::Circle:
                m_out << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << r << "\"/>\n";
                break;
            case Marker::Dot:
                m_out << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << 0.3 * r << "\" fill=\"" << color
                      << "\"/>\n";
                break;
            case Marker::Cross:
                m_out << "<path d=\"M" << x - r << "," << y - r << " L" << x + r << "," << y + r << " M" << x - r
                      << "," << y + r << " L" << x + r << "," << y - r << "\"/>\n";
                break;
            case Marker::Star: {
                m_out << "<polygon points=\"";
                for (int k = 0; k < 10; k++)
                {
                    double angle = -0.5 * 3.14159265358979323846 + k * 3.14159265358979323846 / 5.0;
                    double radius = k % 2 == 0 ? r : 0.4 * r;
                    m_out << (k > 0 ? " " : "") << x + radius * std::cos(angle) << "," << y + radius * std::sin(angle);
                }
                m_out << "\"/>\n";
                break;
            }
            default:
                break;
            }
        }
        m_out << "</g>\n";
    }

    void writeText(const Text &text)
    {
        const char *anchor = text.halign == HAlign::Left ? "start" : text.halign == HAlign::Center ? "middle" : "end";
        double dy = 0.0;
        switch (text.valign)
        {
        case VAlign::Bottom:
            dy = -0.25 * text.fontSize;
            break;
        case VAlign::Center:
            dy = 0.35 * text.fontSize;
            break;
        case VAlign::Top:
            dy = 0.8 * text.fontSize;
            break;
        default:
            break;
        }
        m_out << "<text x=\"" << m_layout.px(text.position.x) << "\" y=\"" << m_layout.py(text.position.y) + dy
              << "\" font-size=\"" << text.fontSize << "\" text-anchor=\"" << anchor << "\" fill=\""
              << svgColor(text.color) << "\" stroke=\"none\">" << escape(text.text) << "</text>\n";
    }

    const Figure &m_fig;
    double m_width;
    double m_height;
    Layout m_layout;
    std::ostringstream m_out;
};

} // namespace

std::string to_svg(const Figure &fig)
{
    return to_svg(fig, fig.width(), fig.height());
}

std::string to_svg(const Figure &fig, double width, double height)
{
    return SvgWriter(fig, width, height, layout(fig, width, height)).write();
}

std::string to_svg(const Figure &fig, double width, double height, const Bounds &view)
{
    return SvgWriter(fig, width, height, layout(fig, width, height, view)).write();
}

} // namespace calfem::vis
