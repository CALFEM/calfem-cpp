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

// Filled polygon in model coordinates (closed). edge is the outline colour (no outline when
// empty); width is the outline width in pixels.
struct Polygon {
    std::vector<Point> points;
    Color fill;
    std::optional<Color> edge;
    double width = 1.0;
    LineStyle style = LineStyle::Solid;
};

using Item = std::variant<Line, Markers, Text, Polygon>;

// ---- Colour maps ----

// matplotlib colour maps: viridis (matplotlib's default), jet (MATLAB's former default) and gray.
enum class Colormap { Viridis, Jet, Gray };

// Colour of cmap at t in [0, 1] (clamped). As matplotlib, the map is a table of 256 colours
// and t selects entry floor(256 t) (entry 255 for t = 1).
Color colormap_color(Colormap cmap, double t);

// Mapping of values to colours of a drawing, shown by colorbar(). A value v gets the colour
// colormap_color(colormap, (v - vmin) / (vmax - vmin)) (0 when vmin == vmax). Filled contours
// set boundaries (the contour levels): the colour bar then shows one band per level interval
// with the colour of its midpoint.
struct ColorScale {
    Colormap colormap = Colormap::Viridis;
    double vmin = 0.0;
    double vmax = 1.0;
    std::vector<double> boundaries;

    // Colour of value (see above).
    Color color(double value) const;
};

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

    // Colour map used by the colour-mapped drawing functions (default viridis).
    Colormap colormap() const;
    void set_colormap(Colormap cmap);

    // Colour scale of the last colour-mapped drawing (matplotlib's current mappable); empty
    // before the first one.
    const std::optional<ColorScale> &color_scale() const;
    void set_color_scale(const ColorScale &scale);

    // Scale of the colour bar at the right of the plot (set by colorbar() to the colour scale
    // at that time, as matplotlib binds a colour bar to the current mappable); empty: no bar.
    const std::optional<ColorScale> &colorbar() const;
    void set_colorbar(const std::optional<ColorScale> &scale);

private:
    int m_number;
    std::string m_title;
    std::vector<Item> m_items;
    std::optional<Bounds> m_axis;
    bool m_equalAspect = true;
    bool m_grid = false;
    double m_width = 800.0;
    double m_height = 600.0;
    Colormap m_colormap = Colormap::Viridis;
    std::optional<ColorScale> m_colorScale;
    std::optional<ColorScale> m_colorbar;
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

// Colour map for the following colour-mapped drawings in the current figure, or fig (as
// matplotlib's set_cmap()).
void colormap(Colormap cmap);
void colormap(Figure &fig, Colormap cmap);

// Shows a colour bar for the last colour-mapped drawing (draw_element_values,
// draw_nodal_values*, eliso2, elval2, ...) in the current figure, or fig, as calfem-python's
// colorbar(). Throws std::runtime_error when the figure has no colour-mapped drawing.
void colorbar();
void colorbar(Figure &fig);

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
// margin; with equal aspect it is widened to give x and y the same scale. With a colour bar
// (Figure::colorbar()) the plot area leaves room for it at the right.

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

    // Colour bar (when colorbar is true): rectangle in pixels, value range from colorbarMin
    // (bottom) to colorbarMax (top) and ticks.
    bool colorbar = false;
    double colorbarLeft = 0.0;
    double colorbarTop = 0.0;
    double colorbarWidth = 0.0;
    double colorbarHeight = 0.0;
    double colorbarMin = 0.0;
    double colorbarMax = 1.0;
    std::vector<double> colorbarTicks;
    double colorbarTickStep = 0.0;

    // Model to pixel coordinates (y pointing down).
    double px(double x) const;
    double py(double y) const;
    // Pixel to model coordinates.
    double mx(double px) const;
    double my(double py) const;
    // Value to pixel y on the colour bar.
    double cy(double value) const;
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

// ---- calfem-python mesh and field drawing functions (draw_*, eliso2, elval2, elflux2) ----
//
// Meshes are given as node coordinates coords (N x 2), element dofs edof (one row per element,
// dofs_per_node dofs per node, numbered node by node: node n has dofs (n-1)*dofs_per_node + 1
// ... n*dofs_per_node, as create_dofs and the mesh module) and the Gmsh element type el_type:
// 2 triangles, 3 quadrilaterals, 16 8-node quadrilaterals (drawn with their corner nodes).
// Colour-mapped drawings use the figure's colour map (colormap()) and set its colour scale
// (colorbar()). Filled regions are Polygon items, contour lines Line items. title sets the
// figure title unless empty.
//
// Deviations from calfem-python:
// - The draw_nodal_values_* functions use edof as node topology (one dof per node) when
//   dofs_per_node is 0, as calfem-python, and convert it with dofs_per_node otherwise
//   (calfem-python uses edof directly, which is only right for one dof per node).
// - Filled contours are drawn as one polygon per triangle and level interval (with an outline
//   in the fill colour, hiding seams between the pieces) instead of matplotlib's merged regions.
// - Gouraud shading (draw_nodal_values_shaded, draw_nodal_values) is approximated by 64 filled
//   contour intervals between the colour limits, coloured with the colour of their midpoint.
// - Arrows (draw_element_flux, elflux2) are polygons in model coordinates with matplotlib's
//   quiver proportions; the shaft width is arrow_width times the width of the mesh (matplotlib:
//   of the axes).
// - draw_displacements uses color (calfem-python always draws (0.3, 0.3, 0.3), now the default)
//   and returns the magnification factor; draw_elements draws show_nodes (ignored in
//   calfem-python). The contour functions return their levels. The axes arguments of
//   calfem-python are replaced by the Figure overloads.

// Mesh lines (filled with face_color if filled); show_nodes marks all nodes in coords.
//
//   draw_mesh(coords, edof, dofs_per_node, el_type, title, color, face_color, filled, show_nodes, node_color)
void draw_mesh(const Matrix &coords, const IntMatrix &edof, int dofs_per_node, int el_type,
               const std::string &title = "", const Color &color = {0.0, 0.0, 0.0},
               const Color &face_color = {0.8, 0.8, 0.8}, bool filled = false, bool show_nodes = false,
               const Color &node_color = {0.0, 0.0, 0.0});
void draw_mesh(Figure &fig, const Matrix &coords, const IntMatrix &edof, int dofs_per_node, int el_type,
               const std::string &title = "", const Color &color = {0.0, 0.0, 0.0},
               const Color &face_color = {0.8, 0.8, 0.8}, bool filled = false, bool show_nodes = false,
               const Color &node_color = {0.0, 0.0, 0.0});

// Elements from element coordinates (one row per element), as closed polygons (closed) or
// polylines.
//
//   draw_elements(ex, ey, title, color, face_color, node_color, line_style, filled, closed, show_nodes)
void draw_elements(const Matrix &ex, const Matrix &ey, const std::string &title = "",
                   const Color &color = {0.0, 0.0, 0.0}, const Color &face_color = {0.8, 0.8, 0.8},
                   const Color &node_color = {0.0, 0.0, 0.0}, LineStyle line_style = LineStyle::Solid,
                   bool filled = false, bool closed = true, bool show_nodes = false);
void draw_elements(Figure &fig, const Matrix &ex, const Matrix &ey, const std::string &title = "",
                   const Color &color = {0.0, 0.0, 0.0}, const Color &face_color = {0.8, 0.8, 0.8},
                   const Color &node_color = {0.0, 0.0, 0.0}, LineStyle line_style = LineStyle::Solid,
                   bool filled = false, bool closed = true, bool show_nodes = false);

// Markers at all element nodes (filled with color if filled).
//
//   draw_node_circles(ex, ey, title, color, face_color, filled, marker_type)
void draw_node_circles(const Matrix &ex, const Matrix &ey, const std::string &title = "",
                       const Color &color = {0.0, 0.0, 0.0}, const Color &face_color = {0.8, 0.8, 0.8},
                       bool filled = false, Marker marker_type = Marker::Circle);
void draw_node_circles(Figure &fig, const Matrix &ex, const Matrix &ey, const std::string &title = "",
                       const Color &color = {0.0, 0.0, 0.0}, const Color &face_color = {0.8, 0.8, 0.8},
                       bool filled = false, Marker marker_type = Marker::Circle);

// Elements filled with the colour of their value (one per element), optionally on the mesh
// displaced by magnfac * displacements (the global displacement vector, dofs_per_node values
// per node, or N x 2 nodal displacements).
//
//   draw_element_values(values, coords, edof, dofs_per_node, el_type, displacements,
//                       draw_elements, draw_undisplaced_mesh, magnfac, title, color, node_color)
void draw_element_values(const ColVec &values, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                         int el_type, const Matrix &displacements = Matrix(), bool draw_elements = true,
                         bool draw_undisplaced_mesh = false, double magnfac = 1.0, const std::string &title = "",
                         const Color &color = {0.0, 0.0, 0.0}, const Color &node_color = {0.0, 0.0, 0.0});
void draw_element_values(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                         int dofs_per_node, int el_type, const Matrix &displacements = Matrix(),
                         bool draw_elements = true, bool draw_undisplaced_mesh = false, double magnfac = 1.0,
                         const std::string &title = "", const Color &color = {0.0, 0.0, 0.0},
                         const Color &node_color = {0.0, 0.0, 0.0});

// Displaced mesh. magnfac < 0 (default) chooses the factor making the largest displacement
// magscale times the model size. Returns the magnification factor.
//
//   magnfac = draw_displacements(a, coords, edof, dofs_per_node, el_type, draw_undisplaced_mesh,
//                                magnfac, magscale, title, color)
double draw_displacements(const Matrix &a, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                          int el_type, bool draw_undisplaced_mesh = false, double magnfac = -1.0,
                          double magscale = 0.1, const std::string &title = "", const Color &color = {0.3, 0.3, 0.3});
double draw_displacements(Figure &fig, const Matrix &a, const Matrix &coords, const IntMatrix &edof,
                          int dofs_per_node, int el_type, bool draw_undisplaced_mesh = false, double magnfac = -1.0,
                          double magscale = 0.1, const std::string &title = "", const Color &color = {0.3, 0.3, 0.3});

// Element flux vectors (flux: one row [qx, qy] per element) as arrows centred at the element
// centroids. scale = arrow length / flux magnitude, by default giving a mean arrow length of 0.8
// times the mean element size. color_by_magnitude colours the arrows by |q| (colour scale).
// Returns the scale.
//
//   scale = draw_element_flux(flux, coords, edof, dofs_per_node, el_type, scale, color,
//                             color_by_magnitude, draw_elements, element_color, arrow_width, title)
double draw_element_flux(const Matrix &flux, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                         int el_type, std::optional<double> scale = std::nullopt, const Color &color = {0.0, 0.0, 0.0},
                         bool color_by_magnitude = false, bool draw_elements = false,
                         const Color &element_color = {0.6, 0.6, 0.6}, double arrow_width = 0.003,
                         const std::string &title = "");
double draw_element_flux(Figure &fig, const Matrix &flux, const Matrix &coords, const IntMatrix &edof,
                         int dofs_per_node, int el_type, std::optional<double> scale = std::nullopt,
                         const Color &color = {0.0, 0.0, 0.0}, bool color_by_magnitude = false,
                         bool draw_elements = false, const Color &element_color = {0.6, 0.6, 0.6},
                         double arrow_width = 0.003, const std::string &title = "");

// Element flux arrows from element coordinates (3- or 4-node elements); plotcolor is the
// colour code. Returns the scale factor (default as draw_element_flux).
//
//   sfac = elflux2(ex, ey, es, plotcolor, sfac)    double sfac = vis::elflux2(ex, ey, es);
double elflux2(const Matrix &ex, const Matrix &ey, const Matrix &es, int plotcolor = 1,
               std::optional<double> sfac = std::nullopt);
double elflux2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &es, int plotcolor = 1,
               std::optional<double> sfac = std::nullopt);

// Filled contours of nodal values (one per node) with levels levels: a number (levels are
// chosen as matplotlib's MaxNLocator(levels + 1)) or the level values (increasing). Regions
// outside the first and last level are not filled. draw_elements draws the mesh on top (needs
// dofs_per_node and el_type). Returns the levels.
//
//   levels = draw_nodal_values_contourf(values, coords, edof, levels, title, dofs_per_node,
//                                       el_type, draw_elements)
std::vector<double> draw_nodal_values_contourf(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                               int levels = 12, const std::string &title = "",
                                               int dofs_per_node = 0, int el_type = 0, bool draw_elements = false);
std::vector<double> draw_nodal_values_contourf(Figure &fig, const ColVec &values, const Matrix &coords,
                                               const IntMatrix &edof, int levels = 12, const std::string &title = "",
                                               int dofs_per_node = 0, int el_type = 0, bool draw_elements = false);
std::vector<double> draw_nodal_values_contourf(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                               const std::vector<double> &levels, const std::string &title = "",
                                               int dofs_per_node = 0, int el_type = 0, bool draw_elements = false);
std::vector<double> draw_nodal_values_contourf(Figure &fig, const ColVec &values, const Matrix &coords,
                                               const IntMatrix &edof, const std::vector<double> &levels,
                                               const std::string &title = "", int dofs_per_node = 0, int el_type = 0,
                                               bool draw_elements = false);

// Contour lines of nodal values, coloured by level (arguments as draw_nodal_values_contourf).
// Returns the levels.
std::vector<double> draw_nodal_values_contour(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                              int levels = 12, const std::string &title = "", int dofs_per_node = 0,
                                              int el_type = 0, bool draw_elements = false);
std::vector<double> draw_nodal_values_contour(Figure &fig, const ColVec &values, const Matrix &coords,
                                              const IntMatrix &edof, int levels = 12, const std::string &title = "",
                                              int dofs_per_node = 0, int el_type = 0, bool draw_elements = false);
std::vector<double> draw_nodal_values_contour(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                                              const std::vector<double> &levels, const std::string &title = "",
                                              int dofs_per_node = 0, int el_type = 0, bool draw_elements = false);
std::vector<double> draw_nodal_values_contour(Figure &fig, const ColVec &values, const Matrix &coords,
                                              const IntMatrix &edof, const std::vector<double> &levels,
                                              const std::string &title = "", int dofs_per_node = 0, int el_type = 0,
                                              bool draw_elements = false);

// Smoothly shaded nodal values (colour limits: the value range).
//
//   draw_nodal_values_shaded(values, coords, edof, title, dofs_per_node, el_type, draw_elements)
void draw_nodal_values_shaded(const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                              const std::string &title = "", int dofs_per_node = 0, int el_type = 0,
                              bool draw_elements = false);
void draw_nodal_values_shaded(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                              const std::string &title = "", int dofs_per_node = 0, int el_type = 0,
                              bool draw_elements = false);

// Smoothly shaded nodal values with the mesh on top (former visvis-based calfem.vis). clim =
// {min, max} sets the colour limits (values outside get the end colours).
//
//   draw_nodal_values(values, coords, edof, dofs_per_node, el_type, clim, draw_elements, title)
void draw_nodal_values(const ColVec &values, const Matrix &coords, const IntMatrix &edof, int dofs_per_node,
                       int el_type, std::optional<std::array<double, 2>> clim = std::nullopt,
                       bool draw_elements = true, const std::string &title = "");
void draw_nodal_values(Figure &fig, const ColVec &values, const Matrix &coords, const IntMatrix &edof,
                       int dofs_per_node, int el_type, std::optional<std::array<double, 2>> clim = std::nullopt,
                       bool draw_elements = true, const std::string &title = "");

// Isolines of element nodal values ed (one row per element, e.g. from extract_eldisp) for
// triangles, quadrilaterals and 8-node elements. isov: number of isolines (levels as
// draw_nodal_values_contour) or their values. plotpar = {linetype, linecolor} draws all lines
// in one colour; by default they are coloured by value. Returns the levels.
//
//   levels = eliso2(ex, ey, ed, isov, plotpar)    vis::eliso2(ex, ey, ed, 10);
std::vector<double> eliso2(const Matrix &ex, const Matrix &ey, const Matrix &ed, int isov = 10,
                           std::optional<std::array<int, 2>> plotpar = std::nullopt);
std::vector<double> eliso2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed, int isov = 10,
                           std::optional<std::array<int, 2>> plotpar = std::nullopt);
std::vector<double> eliso2(const Matrix &ex, const Matrix &ey, const Matrix &ed, const std::vector<double> &isov,
                           std::optional<std::array<int, 2>> plotpar = std::nullopt);
std::vector<double> eliso2(Figure &fig, const Matrix &ex, const Matrix &ey, const Matrix &ed,
                           const std::vector<double> &isov, std::optional<std::array<int, 2>> plotpar = std::nullopt);

// Elements filled with the colour of their value ev (one per element); showMesh draws the
// element edges in black.
//
//   elval2(ex, ey, ev, showMesh)    vis::elval2(ex, ey, ev, true);
void elval2(const Matrix &ex, const Matrix &ey, const ColVec &ev, bool showMesh = false);
void elval2(Figure &fig, const Matrix &ex, const Matrix &ey, const ColVec &ev, bool showMesh = false);

} // namespace calfem::vis
