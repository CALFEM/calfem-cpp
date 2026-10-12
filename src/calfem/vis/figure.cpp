//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/vis.h>

#include "detail/log.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>

namespace calfem::vis {

// ---- Bounds ----

bool Bounds::empty() const
{
    return xmin > xmax || ymin > ymax;
}

void Bounds::include(const Point &p)
{
    if (empty())
    {
        xmin = xmax = p.x;
        ymin = ymax = p.y;
        return;
    }
    xmin = std::min(xmin, p.x);
    xmax = std::max(xmax, p.x);
    ymin = std::min(ymin, p.y);
    ymax = std::max(ymax, p.y);
}

void Bounds::include(const Bounds &b)
{
    if (b.empty())
        return;
    include(Point{b.xmin, b.ymin});
    include(Point{b.xmax, b.ymax});
}

// ---- Figure ----

Figure::Figure(int number) : m_number(number)
{
}

int Figure::number() const
{
    return m_number;
}

const std::string &Figure::title() const
{
    return m_title;
}

void Figure::set_title(const std::string &title)
{
    m_title = title;
}

const std::vector<Item> &Figure::items() const
{
    return m_items;
}

void Figure::add(Item item)
{
    m_items.push_back(std::move(item));
}

void Figure::clear()
{
    m_items.clear();
    m_title.clear();
    m_axis.reset();
    m_colorScale.reset();
    m_colorbar.reset();
}

Bounds Figure::data_bounds() const
{
    Bounds b;
    for (const auto &item : m_items)
    {
        if (auto line = std::get_if<Line>(&item))
            for (const auto &p : line->points)
                b.include(p);
        else if (auto markers = std::get_if<Markers>(&item))
            for (const auto &p : markers->points)
                b.include(p);
        else if (auto text = std::get_if<Text>(&item))
            b.include(text->position);
        else if (auto polygon = std::get_if<Polygon>(&item))
            for (const auto &p : polygon->points)
                b.include(p);
    }
    return b;
}

const std::optional<Bounds> &Figure::axis() const
{
    return m_axis;
}

void Figure::set_axis(const Bounds &limits)
{
    m_axis = limits;
}

void Figure::clear_axis()
{
    m_axis.reset();
}

bool Figure::equal_aspect() const
{
    return m_equalAspect;
}

void Figure::set_equal_aspect(bool equal)
{
    m_equalAspect = equal;
}

bool Figure::grid() const
{
    return m_grid;
}

void Figure::set_grid(bool grid)
{
    m_grid = grid;
}

double Figure::width() const
{
    return m_width;
}

double Figure::height() const
{
    return m_height;
}

void Figure::set_size(double width, double height)
{
    if (width <= 0.0 || height <= 0.0)
        throw std::invalid_argument("Figure::set_size: width and height must be positive.");
    m_width = width;
    m_height = height;
}

Colormap Figure::colormap() const
{
    return m_colormap;
}

void Figure::set_colormap(Colormap cmap)
{
    m_colormap = cmap;
}

const std::optional<ColorScale> &Figure::color_scale() const
{
    return m_colorScale;
}

void Figure::set_color_scale(const ColorScale &scale)
{
    m_colorScale = scale;
}

const std::optional<ColorScale> &Figure::colorbar() const
{
    return m_colorbar;
}

void Figure::set_colorbar(const std::optional<ColorScale> &scale)
{
    m_colorbar = scale;
}

// ---- Figure management ----

namespace {

struct FigureRegistry {
    std::vector<std::unique_ptr<Figure>> figures;
    Figure *current = nullptr;
};

FigureRegistry &registry()
{
    static FigureRegistry r;
    return r;
}

Figure *findFigure(int number)
{
    for (auto &f : registry().figures)
        if (f->number() == number)
            return f.get();
    return nullptr;
}

void removeFigure(Figure *fig)
{
    auto &r = registry();
    r.figures.erase(std::remove_if(r.figures.begin(), r.figures.end(),
                                   [fig](const std::unique_ptr<Figure> &f) { return f.get() == fig; }),
                    r.figures.end());
    if (r.current == fig)
        r.current = r.figures.empty() ? nullptr : r.figures.back().get();
}

} // namespace

Figure &figure()
{
    int number = 1;
    for (auto &f : registry().figures)
        number = std::max(number, f->number() + 1);
    return figure(number);
}

Figure &figure(int number)
{
    auto &r = registry();
    Figure *fig = findFigure(number);
    if (fig == nullptr)
    {
        r.figures.push_back(std::make_unique<Figure>(number));
        fig = r.figures.back().get();
    }
    r.current = fig;
    return *fig;
}

Figure &gcf()
{
    auto &r = registry();
    if (r.current == nullptr)
        return figure();
    return *r.current;
}

void clf()
{
    gcf().clear();
}

void close()
{
    if (registry().current != nullptr)
        removeFigure(registry().current);
}

void close(int number)
{
    if (Figure *fig = findFigure(number))
        removeFigure(fig);
}

void close_all()
{
    registry().figures.clear();
    registry().current = nullptr;
}

std::vector<Figure *> figures()
{
    std::vector<Figure *> result;
    for (auto &f : registry().figures)
        result.push_back(f.get());
    return result;
}

void title(const std::string &title)
{
    gcf().set_title(title);
}

void title(Figure &fig, const std::string &title)
{
    fig.set_title(title);
}

void axis(const std::array<double, 4> &limits)
{
    axis(gcf(), limits);
}

void axis(Figure &fig, const std::array<double, 4> &limits)
{
    if (!(limits[0] < limits[1]) || !(limits[2] < limits[3]))
        throw std::invalid_argument("axis: limits must be {xmin, xmax, ymin, ymax} with xmin < xmax and ymin < ymax.");
    fig.set_axis(Bounds{limits[0], limits[1], limits[2], limits[3]});
}

void axis(const char *option)
{
    axis(gcf(), option);
}

void axis(Figure &fig, const char *option)
{
    std::string o = option != nullptr ? option : "";
    if (o == "equal")
        fig.set_equal_aspect(true);
    else if (o == "auto")
    {
        fig.set_equal_aspect(false);
        fig.clear_axis();
    }
    else if (o == "tight")
        fig.clear_axis();
    else
        throw std::invalid_argument("axis: option must be \"equal\", \"auto\" or \"tight\", got \"" + o + "\".");
}

void colormap(Colormap cmap)
{
    gcf().set_colormap(cmap);
}

void colormap(Figure &fig, Colormap cmap)
{
    fig.set_colormap(cmap);
}

void colorbar()
{
    colorbar(gcf());
}

void colorbar(Figure &fig)
{
    if (!fig.color_scale())
        throw std::runtime_error("colorbar: the figure has no colour-mapped drawing.");
    fig.set_colorbar(fig.color_scale());
}

void show_grid(bool flag)
{
    gcf().set_grid(flag);
}

void show_grid(Figure &fig, bool flag)
{
    fig.set_grid(flag);
}

void text(const std::string &text, const Point &pos, const Color &color, double fontSize)
{
    vis::text(gcf(), text, pos, color, fontSize);
}

void text(Figure &fig, const std::string &text, const Point &pos, const Color &color, double fontSize)
{
    Text t;
    t.position = pos;
    t.text = text;
    t.color = color;
    t.fontSize = fontSize;
    fig.add(t);
}

// ---- Backends ----

namespace {

bool envSet(const char *name)
{
#ifdef _MSC_VER
    char *value = nullptr;
    std::size_t size = 0;
    bool set = _dupenv_s(&value, &size, name) == 0 && value != nullptr;
    std::free(value);
    return set;
#else
    return std::getenv(name) != nullptr;
#endif
}

std::string lowerExtension(const std::string &filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return char(std::tolower(c)); });
    return ext;
}

void writeSvgFile(const std::string &filename, const Figure &fig)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out)
        throw std::runtime_error("save_figure: cannot open " + filename + " for writing.");
    out << to_svg(fig);
    if (!out)
        throw std::runtime_error("save_figure: failed to write " + filename + ".");
}

// Without a backend: write the figures as SVG to the temporary directory and open them with
// the default viewer (no blocking).
void showWithSystemViewer(const std::vector<Figure *> &figs, bool)
{
    std::error_code ec;
    auto dir = std::filesystem::temp_directory_path(ec);
    if (ec)
        dir = std::filesystem::current_path();

    for (auto fig : figs)
    {
        auto path = (dir / ("calfem_figure_" + std::to_string(fig->number()) + ".svg")).string();
        writeSvgFile(path, *fig);
        detail::logMessage("show: figure " + std::to_string(fig->number()) + " written to " + path);

#if defined(_WIN32)
        std::string command = "start \"\" \"" + path + "\"";
#elif defined(__APPLE__)
        std::string command = "open \"" + path + "\"";
#else
        std::string command = "xdg-open \"" + path + "\" >/dev/null 2>&1 &";
#endif
        if (!envSet("CALFEM_NO_VIEWER") && std::system(command.c_str()) != 0)
            detail::logMessage("show: could not open " + path);
    }
}

Backend &currentBackend()
{
    static Backend b{showWithSystemViewer, nullptr};
    return b;
}

} // namespace

void set_backend(const Backend &backend)
{
    currentBackend() = backend;
    if (!currentBackend().show)
        currentBackend().show = showWithSystemViewer;
}

const Backend &backend()
{
    return currentBackend();
}

void save_figure(const std::string &filename)
{
    save_figure(filename, gcf());
}

void save_figure(const std::string &filename, const Figure &fig)
{
    if (lowerExtension(filename) == ".svg")
    {
        writeSvgFile(filename, fig);
        return;
    }

    const auto &b = backend();
    if (b.save && b.save(fig, filename))
        return;

    throw std::runtime_error("save_figure: unsupported file format for " + filename +
                             " (only .svg without a backend such as calfem::vis_qt).");
}

void show()
{
    auto figs = figures();
    if (figs.empty())
        return;
    backend().show(figs, !envSet("CALFEM_NO_BLOCK"));
}

void show_and_wait()
{
    show();
}

// ---- Layout ----

namespace {

// Tick spacing of 1, 2 or 5 times a power of ten giving about maxTicks ticks over range.
double niceStep(double range, double maxTicks)
{
    if (!(range > 0.0) || !(maxTicks >= 1.0))
        return 1.0;
    double raw = range / maxTicks;
    double magnitude = std::pow(10.0, std::floor(std::log10(raw)));
    double residual = raw / magnitude;
    double nice = residual <= 1.0 ? 1.0 : residual <= 2.0 ? 2.0 : residual <= 5.0 ? 5.0 : 10.0;
    return nice * magnitude;
}

std::vector<double> ticks(double lo, double hi, double step)
{
    std::vector<double> result;
    if (!(step > 0.0) || !(hi >= lo))
        return result;
    double first = std::ceil(lo / step - 1e-9) * step;
    for (int i = 0; i < 1000; i++)
    {
        double t = first + i * step;
        if (t > hi + 1e-9 * step)
            break;
        result.push_back(std::abs(t) < 1e-12 * step ? 0.0 : t);
    }
    return result;
}

const double marginLeft = 64.0;
const double marginRight = 20.0;
const double colorbarSpace = 90.0; // colour bar: gap, bar and tick labels
const double colorbarGap = 16.0;
const double colorbarBarWidth = 16.0;
const double marginTop = 36.0;
const double marginBottom = 32.0;

} // namespace

double Layout::px(double x) const
{
    return plotLeft + (x - view.xmin) / (view.xmax - view.xmin) * plotWidth;
}

double Layout::py(double y) const
{
    return plotTop + (view.ymax - y) / (view.ymax - view.ymin) * plotHeight;
}

double Layout::cy(double value) const
{
    return colorbarTop + (colorbarMax - value) / (colorbarMax - colorbarMin) * colorbarHeight;
}

double Layout::mx(double px) const
{
    return view.xmin + (px - plotLeft) / plotWidth * (view.xmax - view.xmin);
}

double Layout::my(double py) const
{
    return view.ymax - (py - plotTop) / plotHeight * (view.ymax - view.ymin);
}

Layout layout(const Figure &fig, double width, double height)
{
    Bounds view;
    if (fig.axis())
        view = *fig.axis();
    else
    {
        view = fig.data_bounds();
        if (view.empty())
            view = Bounds{0.0, 1.0, 0.0, 1.0};

        // Degenerate extents (a single point or a horizontal/vertical line)
        double size = std::max(view.xmax - view.xmin, view.ymax - view.ymin);
        if (!(size > 0.0))
            size = std::max(1.0, std::abs(view.xmin) + std::abs(view.ymin));
        if (!(view.xmax > view.xmin))
        {
            view.xmin -= 0.5 * size;
            view.xmax += 0.5 * size;
        }
        if (!(view.ymax > view.ymin))
        {
            view.ymin -= 0.5 * size;
            view.ymax += 0.5 * size;
        }

        // 5% margin, as matplotlib's default axes margins
        double mx = 0.05 * (view.xmax - view.xmin), my = 0.05 * (view.ymax - view.ymin);
        view = Bounds{view.xmin - mx, view.xmax + mx, view.ymin - my, view.ymax + my};
    }
    return layout(fig, width, height, view);
}

Layout layout(const Figure &fig, double width, double height, const Bounds &view)
{
    Layout l;
    l.colorbar = fig.colorbar().has_value();
    l.plotLeft = marginLeft;
    l.plotTop = marginTop;
    l.plotWidth = std::max(10.0, width - marginLeft - marginRight - (l.colorbar ? colorbarSpace : 0.0));
    l.plotHeight = std::max(10.0, height - marginTop - marginBottom);
    l.view = view;

    if (l.colorbar)
    {
        const ColorScale &scale = *fig.colorbar();
        l.colorbarLeft = l.plotLeft + l.plotWidth + colorbarGap;
        l.colorbarTop = l.plotTop;
        l.colorbarWidth = colorbarBarWidth;
        l.colorbarHeight = l.plotHeight;

        if (scale.boundaries.size() >= 2)
        {
            // Filled contours: ticks at the levels (every k-th level when there are many)
            l.colorbarMin = scale.boundaries.front();
            l.colorbarMax = scale.boundaries.back();
            std::size_t n = scale.boundaries.size();
            std::size_t every = std::max<std::size_t>(1, (n + 11) / 12);
            for (std::size_t i = 0; i < n; i += every)
                l.colorbarTicks.push_back(scale.boundaries[i]);

            // Label precision: the decimals the levels need (0.25 must not become 0.2), relative
            // to the level range
            double range = std::max(std::abs(l.colorbarMax - l.colorbarMin), 1e-300);
            double step = std::pow(10.0, std::floor(std::log10(range)));
            for (int k = 0; k < 12; k++)
            {
                bool exact = true;
                for (double t : scale.boundaries)
                    exact = exact && std::abs(t / step - std::round(t / step)) < 1e-6;
                if (exact)
                    break;
                step /= 10.0;
            }
            l.colorbarTickStep = step;
        }
        else
        {
            l.colorbarMin = scale.vmin;
            l.colorbarMax = scale.vmax;
            l.colorbarTickStep = niceStep(l.colorbarMax - l.colorbarMin, l.colorbarHeight / 40.0);
            l.colorbarTicks = ticks(l.colorbarMin, l.colorbarMax, l.colorbarTickStep);
        }

        // A constant field: show a unit range around the value
        if (!(l.colorbarMax > l.colorbarMin))
        {
            l.colorbarMin -= 0.5;
            l.colorbarMax += 0.5;
            l.colorbarTickStep = niceStep(1.0, l.colorbarHeight / 40.0);
            l.colorbarTicks = ticks(l.colorbarMin, l.colorbarMax, l.colorbarTickStep);
        }
    }

    if (fig.equal_aspect())
    {
        // Widen the smaller extent so that one model unit has the same length in x and y
        double scale = std::max((view.xmax - view.xmin) / l.plotWidth, (view.ymax - view.ymin) / l.plotHeight);
        double cx = 0.5 * (view.xmin + view.xmax), cy = 0.5 * (view.ymin + view.ymax);
        double hw = 0.5 * scale * l.plotWidth, hh = 0.5 * scale * l.plotHeight;
        l.view = Bounds{cx - hw, cx + hw, cy - hh, cy + hh};
    }

    l.tickStepX = niceStep(l.view.xmax - l.view.xmin, l.plotWidth / 90.0);
    l.tickStepY = niceStep(l.view.ymax - l.view.ymin, l.plotHeight / 60.0);
    l.xticks = ticks(l.view.xmin, l.view.xmax, l.tickStepX);
    l.yticks = ticks(l.view.ymin, l.view.ymax, l.tickStepY);

    return l;
}

std::string format_tick(double value, double step)
{
    if (std::abs(value) < 1e-12 * std::abs(step))
        value = 0.0;

    double magnitude = std::max(std::abs(value), std::abs(step));
    char buffer[64];
    if (magnitude != 0.0 && (magnitude >= 1e5 || magnitude < 1e-3))
    {
        // Scientific notation with enough digits to separate neighbouring ticks
        int digits = std::max(0, int(std::floor(std::log10(magnitude))) - int(std::floor(std::log10(step))));
        std::snprintf(buffer, sizeof(buffer), "%.*e", std::min(digits, 15), value);
    }
    else
    {
        int decimals = step > 0.0 ? std::max(0, -int(std::floor(std::log10(step) + 1e-9))) : 0;
        std::snprintf(buffer, sizeof(buffer), "%.*f", std::min(decimals, 15), value);
    }
    return buffer;
}

} // namespace calfem::vis
