//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Visualisation (port of calfem-python's calfem.vis_mpl).
//
// The drawing functions add lines, markers and text to a Figure. A figure is written to SVG
// with save_figure() or shown with show(). show() uses the registered backend: include
// <calfem/vis_qt.h> and link calfem::vis_qt for interactive Qt windows; without a backend the
// figures are written as SVG files to the temporary directory and opened with the default viewer.
//
//   vis::figure(1);
//   vis::eldraw2(ex, ey, {1, 2, 1});
//   double sfac = vis::eldisp2(ex, ey, ed, {2, 1, 1});
//   vis::title("Displacements");
//   vis::save_figure("frame.svg");
//   vis::show();
//
// As in matplotlib, the functions without a Figure argument draw in the current figure
// (figure(), gcf()). The current-figure functions are not thread-safe; every function also has
// an overload taking the Figure to draw in.

#pragma once

#include <calfem/types.h>

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace calfem::vis {

// ---- Figure model ----

struct Point {
    double x = 0.0;
    double y = 0.0;
};

// RGB colour, components in [0, 1] (as matplotlib colour tuples).
struct Color {
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
};

enum class LineStyle { Solid, Dashed, Dotted };

enum class Marker { None, Circle, Cross, Star, Dot };

enum class HAlign { Left, Center, Right };

enum class VAlign { Baseline, Bottom, Center, Top };

// Polyline in model coordinates; closed adds the segment from the last to the first point.
// width is in pixels (points at 72 dpi).
struct Line {
    std::vector<Point> points;
    Color color;
    LineStyle style = LineStyle::Solid;
    double width = 1.0;
    bool closed = false;
};

// Markers at points; size is the marker diameter in pixels.
struct Markers {
    std::vector<Point> points;
    Marker marker = Marker::Circle;
    Color color;
    double size = 6.0;
    bool filled = false;
};

// Text anchored at position (model coordinates); fontSize in points.
struct Text {
    Point position;
    std::string text;
    Color color;
    double fontSize = 10.0;
    HAlign halign = HAlign::Left;
    VAlign valign = VAlign::Baseline;
};

using Item = std::variant<Line, Markers, Text>;

// Axis-aligned rectangle in model coordinates. Default constructed bounds are empty.
struct Bounds {
    double xmin = 0.0;
    double xmax = -1.0;
    double ymin = 0.0;
    double ymax = -1.0;

    bool empty() const;
    void include(const Point &p);
    void include(const Bounds &b);
};

class Figure {
public:
    explicit Figure(int number = 1);

    int number() const;

    const std::string &title() const;
    void set_title(const std::string &title);

    const std::vector<Item> &items() const;
    void add(Item item);

    // Removes all items, the title and the axis limits.
    void clear();

    // Bounds of all item points (empty when there are no items).
    Bounds data_bounds() const;

    // Fixed axis limits; without them the view fits the data with a 5% margin.
    const std::optional<Bounds> &axis() const;
    void set_axis(const Bounds &limits);
    void clear_axis();

    // Equal scaling of x and y (default on, as calfem-python's drawing functions set it).
    bool equal_aspect() const;
    void set_equal_aspect(bool equal);

    bool grid() const;
    void set_grid(bool grid);

    // Size in pixels used by save_figure() and as the initial window size.
    double width() const;
    double height() const;
    void set_size(double width, double height);

private:
    int m_number;
    std::string m_title;
    std::vector<Item> m_items;
    std::optional<Bounds> m_axis;
    bool m_equalAspect = true;
    bool m_grid = false;
    double m_width = 800.0;
    double m_height = 600.0;
};

// ---- Figure management (current figure, as matplotlib.pyplot) ----

// figure(): creates a new figure and makes it current.
// figure(number): makes figure number current, creating it if needed.
Figure &figure();
Figure &figure(int number);

// Current figure (created if there is none).
Figure &gcf();

// Clears the current figure.
void clf();

// Closes (removes) the current figure, figure number, or all figures.
void close();
void close(int number);
void close_all();

// All open figures in creation order.
std::vector<Figure *> figures();

void title(const std::string &title);
void title(Figure &fig, const std::string &title);

// axis({xmin, xmax, ymin, ymax}) fixes the axis limits; axis("equal") sets equal scaling,
// axis("auto") removes fixed limits and equal scaling, axis("tight") fits the data.
void axis(const std::array<double, 4> &limits);
void axis(Figure &fig, const std::array<double, 4> &limits);
void axis(const char *option);
void axis(Figure &fig, const char *option);

void show_grid(bool flag = true);
void show_grid(Figure &fig, bool flag);

// Text at pos (model coordinates), as calfem-python's text()/add_text().
void text(const std::string &text, const Point &pos, const Color &color = {}, double fontSize = 10.0);
void text(Figure &fig, const std::string &text, const Point &pos, const Color &color = {}, double fontSize = 10.0);

// ---- Output ----

// SVG document of the figure (width x height pixels, defaults to the figure size), showing the
// default view or the model region view (see layout()).
std::string to_svg(const Figure &fig);
std::string to_svg(const Figure &fig, double width, double height);
std::string to_svg(const Figure &fig, double width, double height, const Bounds &view);

// Writes the current figure, or fig, to filename. The format follows the extension: .svg is
// always available, other formats (.png, .jpg, .pdf) need a backend (calfem::vis_qt).
// Throws std::runtime_error when the file cannot be written.
void save_figure(const std::string &filename);
void save_figure(const std::string &filename, const Figure &fig);

// Shows all open figures. With the Qt backend show() opens one window per figure and blocks
// until they are closed (unless the environment variable CALFEM_NO_BLOCK is set);
// show_and_wait() is an alias (calfem-python name).
void show();
void show_and_wait();

// ---- Backends ----
//
// A backend shows figures and saves formats other than SVG. calfem::vis_qt installs itself
// when <calfem/vis_qt.h> is included; set_backend() replaces the current backend.

struct Backend {
    // Shows the figures; block: wait until the windows are closed.
    std::function<void(const std::vector<Figure *> &figures, bool block)> show;
    // Saves fig to filename; returns false if the format is not supported.
    std::function<bool(const Figure &fig, const std::string &filename)> save;
};

void set_backend(const Backend &backend);
const Backend &backend();

// ---- Layout (used by renderers) ----
//
// Plot area, visible model region and axis ticks for drawing fig in a width x height pixel
// area. The visible region is view if given, else the axis limits or the data bounds with a 5%
// margin; with equal aspect it is widened to give x and y the same scale.

struct Layout {
    double plotLeft = 0.0;
    double plotTop = 0.0;
    double plotWidth = 0.0;
    double plotHeight = 0.0;
    Bounds view;
    std::vector<double> xticks;
    std::vector<double> yticks;
    double tickStepX = 0.0;
    double tickStepY = 0.0;

    // Model to pixel coordinates (y pointing down).
    double px(double x) const;
    double py(double y) const;
    // Pixel to model coordinates.
    double mx(double px) const;
    double my(double py) const;
};

Layout layout(const Figure &fig, double width, double height);
Layout layout(const Figure &fig, double width, double height, const Bounds &view);

// Tick label for value on an axis with tick spacing step.
std::string format_tick(double value, double step);

// ---- calfem-python drawing functions (2D, element coordinates) ----
//
// ex, ey: element node coordinates, one row per element. plotpar codes as in CALFEM:
//   linetype   1 solid, 2 dashed, 3 dotted
//   linecolor  1 black, 2 blue, 3 magenta, 4 red
//   nodemark   0 none, 1 circle, 2 cross (eldraw2) or star (eldisp2, dispbeam2), 3 dot
// Invalid codes throw std::invalid_argument.

// Undeformed mesh of bars/beams (2 nodes), triangles (3), quadrilaterals (4) or 8-node
// elements; elnum (optional, one per element) is written at the element centres.
//
//   eldraw2(ex, ey, plotpar, elnum)    vis::eldraw2(ex, ey, {1, 2, 1}, elnum);
void eldraw2(const Matrix &ex, const Matrix &ey, const std::array<int, 3> &plotpar = {1, 2, 1},
             const IntRowVec &elnum = IntRowVec());
void eldraw2(Figure &fig, const Matrix &ex, const Matrix &ey, const std::array<int, 3> &plotpar = {1, 2, 1},
             const IntRowVec &elnum = IntRowVec());

// Deformed mesh, displacements scaled by sfac (default 0.1 times the largest model dimension
// divided by the largest displacement). Supports bars (ed: 4 columns), beams (6 columns, drawn
// with beam2crd), triangles, quadrilaterals and 8-node elements. Returns sfac.
//
//   sfac = eldisp2(ex, ey, ed, plotpar, sfac)    double sfac = vis::eldisp2(ex, ey, ed, {2, 1, 1});
double eldisp2(const Matrix &ex, const Matrix &ey, const Matrix &ed, const std::array<int, 3> &plotpar = {2, 1, 1},
               std::optional<double> sfac = std::nullopt);
double eldisp2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed,
               const std::array<int, 3> &plotpar = {2, 1, 1}, std::optional<double> sfac = std::nullopt);

// Displacement diagram of one beam element from the local displacements edi (one row
// [u, v] per evaluation point, e.g. from beam2s), scaled by sfac (default 0.1 L / max |edi|).
// Returns sfac.
//
// Deviation from calfem-python: the line type of plotpar is used (calfem-python always draws
// a solid line); a zero displacement field gives sfac = 1 instead of a division by zero.
//
//   sfac = dispbeam2(ex, ey, edi, plotpar, sfac)    double sfac = vis::dispbeam2(ex, ey, edi);
double dispbeam2(const RowVec &ex, const RowVec &ey, const Matrix &edi, const std::array<int, 3> &plotpar = {2, 1, 1},
                 std::optional<double> sfac = std::nullopt);
double dispbeam2(Figure &fig, const RowVec &ex, const RowVec &ey, const Matrix &edi,
                 const std::array<int, 3> &plotpar = {2, 1, 1}, std::optional<double> sfac = std::nullopt);

// Section force diagram of one bar or beam element. es holds one section force per evaluation
// point (a vector), eci their local x-coordinates (default: evenly spaced over the element).
// plotpar = {diagram colour, element colour}; sfac defaults to 0.2 L / max |es|. Returns sfac.
//
// Deviation from calfem-python: a zero section force gives sfac = 1 instead of a division by zero.
//
//   sfac = secforce2(ex, ey, es, plotpar, sfac, eci)    double sfac = vis::secforce2(ex, ey, es.col(2));
double secforce2(const RowVec &ex, const RowVec &ey, const Matrix &es, const std::array<int, 2> &plotpar = {2, 1},
                 std::optional<double> sfac = std::nullopt, const Matrix &eci = Matrix());
double secforce2(Figure &fig, const RowVec &ex, const RowVec &ey, const Matrix &es,
                 const std::array<int, 2> &plotpar = {2, 1}, std::optional<double> sfac = std::nullopt,
                 const Matrix &eci = Matrix());

// Scale bar of length magnitude * sfac labelled with magnitude. magnitude is {N} (drawn at
// (0, -0.5)) or {N, x, y}; plotpar is the colour code.
//
//   scalgraph2(sfac, magnitude, plotpar)    vis::scalgraph2(sfac, {3e4, 1.5, 0.0});
void scalgraph2(double sfac, const std::vector<double> &magnitude, int plotpar = 2);
void scalgraph2(Figure &fig, double sfac, const std::vector<double> &magnitude, int plotpar = 2);

} // namespace calfem::vis
