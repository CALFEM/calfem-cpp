//
// Tests of the calfem::vis figure model, figure management, layout and SVG output.
//

#include <calfem/vis.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace calfem;
using Catch::Approx;

namespace {

int count(const std::string &text, const std::string &pattern)
{
    int n = 0;
    for (auto pos = text.find(pattern); pos != std::string::npos; pos = text.find(pattern, pos + 1))
        n++;
    return n;
}

Matrix frameEx()
{
    Matrix ex(3, 2);
    ex << 0.0, 0.0, 0.0, 6.0, 6.0, 6.0;
    return ex;
}

Matrix frameEy()
{
    Matrix ey(3, 2);
    ey << 0.0, 4.0, 4.0, 4.0, 4.0, 0.0;
    return ey;
}

} // namespace

TEST_CASE("Figure management follows matplotlib.pyplot", "[vis]")
{
    vis::close_all();
    REQUIRE(vis::figures().empty());

    vis::Figure &f1 = vis::figure(1);
    vis::Figure &f3 = vis::figure(3);
    CHECK(&vis::gcf() == &f3);

    CHECK(&vis::figure(1) == &f1); // selects the existing figure
    CHECK(&vis::gcf() == &f1);

    vis::Figure &f4 = vis::figure(); // next free number
    CHECK(f4.number() == 4);
    CHECK(vis::figures().size() == 3);

    vis::title("Title");
    CHECK(f4.title() == "Title");

    vis::close(); // closes the current figure (4)
    CHECK(vis::figures().size() == 2);
    CHECK(vis::gcf().number() == 3);

    vis::close(1);
    REQUIRE(vis::figures().size() == 1);
    CHECK(vis::figures()[0]->number() == 3);

    vis::close_all();
    CHECK(vis::figures().empty());
    CHECK(vis::gcf().number() == 1); // gcf creates a figure
    vis::close_all();
}

TEST_CASE("Drawing in the current figure", "[vis]")
{
    vis::close_all();
    vis::figure(1);
    vis::eldraw2(frameEx(), frameEy());
    vis::text("label", {1.0, 2.0});

    const auto &items = vis::gcf().items();
    REQUIRE(items.size() == 5); // 3 element lines, node markers, text
    CHECK(std::holds_alternative<vis::Markers>(items[3]));
    CHECK(std::get<vis::Text>(items[4]).text == "label");

    vis::clf();
    CHECK(vis::gcf().items().empty());
    vis::close_all();
}

TEST_CASE("Axis limits and options", "[vis]")
{
    vis::Figure fig;
    vis::eldraw2(fig, frameEx(), frameEy());

    vis::Bounds b = fig.data_bounds();
    CHECK(b.xmin == 0.0);
    CHECK(b.xmax == 6.0);
    CHECK(b.ymin == 0.0);
    CHECK(b.ymax == 4.0);

    vis::axis(fig, {-1.5, 7.5, -0.5, 5.5});
    REQUIRE(fig.axis().has_value());
    CHECK(fig.axis()->xmin == -1.5);
    CHECK(fig.axis()->ymax == 5.5);

    vis::axis(fig, "auto");
    CHECK_FALSE(fig.axis().has_value());
    CHECK_FALSE(fig.equal_aspect());
    vis::axis(fig, "equal");
    CHECK(fig.equal_aspect());

    CHECK_THROWS_AS(vis::axis(fig, {1.0, 0.0, 0.0, 1.0}), std::invalid_argument);
    CHECK_THROWS_AS(vis::axis(fig, "square"), std::invalid_argument);
}

TEST_CASE("Layout keeps equal scaling and maps model to pixels", "[vis]")
{
    vis::Figure fig;
    vis::eldraw2(fig, frameEx(), frameEy());

    vis::Layout l = vis::layout(fig, 800.0, 600.0);

    // Equal aspect: same model units per pixel in x and y
    double sx = (l.view.xmax - l.view.xmin) / l.plotWidth;
    double sy = (l.view.ymax - l.view.ymin) / l.plotHeight;
    CHECK(sx == Approx(sy));

    // The data (with margin) is inside the view
    CHECK(l.view.xmin < 0.0);
    CHECK(l.view.xmax > 6.0);
    CHECK(l.view.ymin < 0.0);
    CHECK(l.view.ymax > 4.0);

    // y axis points up in model coordinates, down in pixels
    CHECK(l.py(4.0) < l.py(0.0));
    CHECK(l.mx(l.px(2.5)) == Approx(2.5));
    CHECK(l.my(l.py(-1.25)) == Approx(-1.25));

    // Ticks are inside the view and evenly spaced
    REQUIRE(l.xticks.size() >= 2);
    for (double t : l.xticks)
    {
        CHECK(t >= l.view.xmin - 1e-12);
        CHECK(t <= l.view.xmax + 1e-12);
    }
    CHECK(l.xticks[1] - l.xticks[0] == Approx(l.tickStepX));

    // Explicit view, no equal aspect
    fig.set_equal_aspect(false);
    vis::Layout l2 = vis::layout(fig, 800.0, 600.0, vis::Bounds{0.0, 10.0, 0.0, 1.0});
    CHECK(l2.view.xmax == 10.0);
    CHECK(l2.view.ymax == 1.0);
    CHECK(l2.px(0.0) == Approx(l2.plotLeft));
    CHECK(l2.py(0.0) == Approx(l2.plotTop + l2.plotHeight));
}

TEST_CASE("Layout of an empty figure and of a single point", "[vis]")
{
    vis::Figure empty;
    vis::Layout l = vis::layout(empty, 400.0, 300.0);
    CHECK(l.view.xmax > l.view.xmin);
    CHECK(l.view.ymax > l.view.ymin);

    vis::Figure point;
    vis::text(point, "x", {2.0, 3.0});
    vis::Layout lp = vis::layout(point, 400.0, 300.0);
    CHECK(lp.view.xmin < 2.0);
    CHECK(lp.view.xmax > 2.0);
    CHECK(lp.view.ymin < 3.0);
    CHECK(lp.view.ymax > 3.0);
}

TEST_CASE("Tick labels", "[vis]")
{
    CHECK(vis::format_tick(0.0, 1.0) == "0");
    CHECK(vis::format_tick(2.0, 1.0) == "2");
    CHECK(vis::format_tick(1.5, 0.5) == "1.5");
    CHECK(vis::format_tick(-0.25, 0.05) == "-0.25");
    CHECK(vis::format_tick(3e-17, 0.1) == "0.0");
    CHECK(vis::format_tick(2e5, 1e5) == "2e+05");
    CHECK(vis::format_tick(1.5e-4, 5e-5) == "1.5e-04");
}

TEST_CASE("SVG output", "[vis]")
{
    vis::Figure fig;
    vis::eldraw2(fig, frameEx(), frameEy(), {2, 4, 1}, IntRowVec::LinSpaced(3, 1, 3));
    vis::title(fig, "Frame <1> & \"2\"");
    fig.set_grid(true);

    std::string svg = vis::to_svg(fig, 640.0, 480.0);

    CHECK(svg.rfind("<?xml", 0) == 0);
    CHECK(svg.find("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"480\"") != std::string::npos);
    CHECK(svg.find("</svg>") != std::string::npos);
    CHECK(count(svg, "<polyline") == 3);              // one per element
    CHECK(count(svg, "stroke=\"#ff0000\"") >= 3);      // red lines (plotpar[1] = 4)
    CHECK(count(svg, "stroke-dasharray") == 3);         // dashed (plotpar[0] = 2)
    CHECK(count(svg, "<circle") == 6);                 // node markers
    CHECK(svg.find("Frame &lt;1&gt; &amp; &quot;2&quot;") != std::string::npos); // escaped title
    CHECK(svg.find(">2</text>") != std::string::npos); // element number

    // Numbers are written with '.' as decimal separator whatever the global locale
    CHECK(svg.find("NaN") == std::string::npos);
    CHECK(svg.find("inf") == std::string::npos);
}

TEST_CASE("save_figure writes SVG and needs a backend for other formats", "[vis]")
{
    vis::Figure fig;
    vis::eldraw2(fig, frameEx(), frameEy());

    auto path = std::filesystem::temp_directory_path() / "calfem_test_figure.svg";
    vis::save_figure(path.string(), fig);

    std::ifstream file(path);
    REQUIRE(file.good());
    std::stringstream content;
    content << file.rdbuf();
    CHECK(content.str() == vis::to_svg(fig));
    file.close();
    std::filesystem::remove(path);

    vis::Backend saved = vis::backend();
    vis::set_backend(vis::Backend{nullptr, nullptr});
    CHECK_THROWS_AS(vis::save_figure("figure.png", fig), std::runtime_error);

    // A backend handles other formats
    std::string savedName;
    vis::set_backend(vis::Backend{nullptr, [&](const vis::Figure &, const std::string &filename) {
                                      savedName = filename;
                                      return true;
                                  }});
    vis::save_figure("figure.png", fig);
    CHECK(savedName == "figure.png");

    vis::set_backend(saved);
}

TEST_CASE("show passes the open figures to the backend", "[vis]")
{
    vis::Backend saved = vis::backend();

    std::vector<int> shown;
    vis::set_backend(vis::Backend{[&](const std::vector<vis::Figure *> &figs, bool) {
                                      for (auto f : figs)
                                          shown.push_back(f->number());
                                  },
                                  nullptr});

    vis::close_all();
    vis::show(); // no figures: backend not called
    CHECK(shown.empty());

    vis::figure(2);
    vis::figure(5);
    vis::show_and_wait();
    CHECK(shown == std::vector<int>{2, 5});

    vis::close_all();
    vis::set_backend(saved);
}

TEST_CASE("Invalid plotpar codes and inputs throw", "[vis]")
{
    vis::Figure fig;
    Matrix ex = frameEx(), ey = frameEy();

    CHECK_THROWS_AS(vis::eldraw2(fig, ex, ey, {4, 1, 1}), std::invalid_argument);
    CHECK_THROWS_AS(vis::eldraw2(fig, ex, ey, {1, 5, 1}), std::invalid_argument);
    CHECK_THROWS_AS(vis::eldraw2(fig, ex, ey, {1, 1, 3}), std::invalid_argument);
    CHECK_THROWS_AS(vis::eldraw2(fig, ex, ey, {1, 1, 1}, IntRowVec::Ones(2)), std::invalid_argument);
    CHECK_THROWS_AS(vis::eldraw2(fig, ex, Matrix(ey.leftCols(1))), std::invalid_argument);

    CHECK_THROWS_AS(vis::eldisp2(fig, ex, ey, Matrix::Zero(3, 5)), std::invalid_argument);
    CHECK_THROWS_AS(vis::eldisp2(fig, ex, ey, Matrix::Zero(2, 6)), std::invalid_argument);

    RowVec bx(2), by(2);
    bx << 0.0, 3.0;
    by << 0.0, 4.0;
    CHECK_THROWS_AS(vis::dispbeam2(fig, bx, by, Matrix::Zero(5, 3)), std::invalid_argument);
    CHECK_THROWS_AS(vis::secforce2(fig, bx, by, Matrix::Zero(3, 2)), std::invalid_argument);
    CHECK_THROWS_AS(vis::secforce2(fig, bx, by, ColVec::Ones(3), {2, 1}, 1.0, ColVec::Ones(4)),
                    std::invalid_argument);
    CHECK_THROWS_AS(vis::scalgraph2(fig, 1.0, {1.0, 2.0}), std::invalid_argument);
    CHECK_THROWS_AS(vis::scalgraph2(fig, 1.0, {1.0}, 7), std::invalid_argument);
}

TEST_CASE("Zero displacements and section forces give a unit scale factor", "[vis]")
{
    vis::Figure fig;
    RowVec bx(2), by(2);
    bx << 0.0, 3.0;
    by << 0.0, 4.0;

    CHECK(vis::eldisp2(fig, frameEx(), frameEy(), Matrix::Zero(3, 6)) == 1.0);
    CHECK(vis::dispbeam2(fig, bx, by, Matrix::Zero(5, 2)) == 1.0);
    CHECK(vis::secforce2(fig, bx, by, ColVec::Zero(5)) == 1.0);
}

TEST_CASE("scalgraph2 with a single magnitude uses the default position", "[vis]")
{
    vis::Figure fig;
    vis::scalgraph2(fig, 2.0, {0.5});

    REQUIRE(fig.items().size() == 4);
    const auto &line = std::get<vis::Line>(fig.items()[0]);
    CHECK(line.points[0].x == 0.0);
    CHECK(line.points[0].y == -0.5);
    CHECK(line.points[1].x == 1.0);
    CHECK(std::get<vis::Text>(fig.items()[3]).text == "0.5");
}
