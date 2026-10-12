//
// Tests of the Qt backend for calfem::vis (calfem::vis_qt). Run with QT_QPA_PLATFORM=minimal
// (set by tests/CMakeLists.txt), so no windows appear.
//

#include <calfem/vis_qt.h>

#include "figure_widget.h"
#include "painter.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <QApplication>
#include <QFileInfo>
#include <QImage>
#include <QPixmap>
#include <QWidget>

#include <filesystem>
#include <memory>

using namespace calfem;
using Catch::Approx;

namespace {

vis::Figure frameFigure()
{
    Matrix ex(3, 2), ey(3, 2);
    ex << 0.0, 0.0, 0.0, 6.0, 6.0, 6.0;
    ey << 0.0, 4.0, 4.0, 4.0, 4.0, 0.0;

    vis::Figure fig(1);
    vis::eldraw2(fig, ex, ey, {1, 2, 1}, IntRowVec::LinSpaced(3, 1, 3));
    vis::title(fig, "Frame");
    return fig;
}

// Number of pixels close to colour (r, g, b)
int countPixels(const QImage &image, int r, int g, int b)
{
    int n = 0;
    QImage img = image.convertToFormat(QImage::Format_RGB32);
    for (int y = 0; y < img.height(); y++)
        for (int x = 0; x < img.width(); x++)
        {
            QColor c = img.pixelColor(x, y);
            if (std::abs(c.red() - r) < 60 && std::abs(c.green() - g) < 60 && std::abs(c.blue() - b) < 60)
                n++;
        }
    return n;
}

} // namespace

TEST_CASE("Including vis_qt.h installs the Qt backend", "[vis][qt]")
{
    CHECK(vis::backend().show);
    CHECK(vis::backend().save);
}

TEST_CASE("figure_widget draws the figure", "[vis][qt]")
{
    vis::Figure fig = frameFigure();

    std::unique_ptr<QWidget> widget(vis::figure_widget(fig));
    REQUIRE(QApplication::instance() != nullptr);
    widget->resize(400, 300);

    QImage image = widget->grab().toImage();
    CHECK(image.width() >= 400);
    CHECK(countPixels(image, 0, 0, 255) > 100); // blue element lines (plotpar[1] = 2)
    CHECK(countPixels(image, 255, 255, 255) > image.width() * image.height() / 2);
}

TEST_CASE("FigureWidget zooms and fits the view", "[vis][qt]")
{
    vis::Figure fig = frameFigure();
    vis::qt::ensureApplication(); // ctest runs each case in its own process
    vis::qt::FigureWidget widget(fig);
    widget.resize(400, 300);

    vis::Layout l0 = widget.currentLayout();
    double w0 = l0.view.xmax - l0.view.xmin;

    widget.zoom(2.0);
    vis::Layout l1 = widget.currentLayout();
    CHECK(l1.view.xmax - l1.view.xmin == Approx(0.5 * w0));
    // Zoom around the centre keeps the centre
    CHECK(0.5 * (l1.view.xmin + l1.view.xmax) == Approx(0.5 * (l0.view.xmin + l0.view.xmax)));

    // Zoom around a point keeps that point under the cursor
    QPointF p(l1.plotLeft + 50.0, l1.plotTop + 40.0);
    double mx = l1.mx(p.x()), my = l1.my(p.y());
    widget.zoom(1.5, p);
    vis::Layout l2 = widget.currentLayout();
    CHECK(l2.mx(p.x()) == Approx(mx));
    CHECK(l2.my(p.y()) == Approx(my));

    widget.fit();
    vis::Layout l3 = widget.currentLayout();
    CHECK(l3.view.xmin == Approx(l0.view.xmin));
    CHECK(l3.view.xmax == Approx(l0.view.xmax));
}

TEST_CASE("save_figure writes PNG and PDF through the Qt backend", "[vis][qt]")
{
    vis::Figure fig = frameFigure();
    fig.set_size(320.0, 240.0);

    auto dir = std::filesystem::temp_directory_path();
    auto png = (dir / "calfem_test_figure.png").string();
    auto pdf = (dir / "calfem_test_figure.pdf").string();

    vis::save_figure(png, fig);
    QImage image(QString::fromStdString(png));
    REQUIRE_FALSE(image.isNull());
    CHECK(image.width() == 640); // rendered at twice the figure size
    CHECK(image.height() == 480);
    CHECK(countPixels(image, 0, 0, 255) > 100);

    vis::save_figure(pdf, fig);
    CHECK(QFileInfo(QString::fromStdString(pdf)).size() > 0);

    CHECK_FALSE(vis::qt::saveFigure(fig, "figure.xyz"));

    std::filesystem::remove(png);
    std::filesystem::remove(pdf);
}

TEST_CASE("show opens one window per figure", "[vis][qt]")
{
    vis::close_all();
    vis::figure(1);
    vis::title("first");
    vis::figure(2);

    int before = 0;
    for (auto w : QApplication::topLevelWidgets())
        if (dynamic_cast<vis::qt::FigureWindow *>(w))
            before++;

    vis::backend().show(vis::figures(), false); // without blocking

    std::vector<vis::qt::FigureWindow *> windows;
    for (auto w : QApplication::topLevelWidgets())
        if (auto fw = dynamic_cast<vis::qt::FigureWindow *>(w))
            windows.push_back(fw);

    CHECK(int(windows.size()) - before == 2);
    bool found = false;
    for (auto w : windows)
        found = found || w->windowTitle() == "Figure 1: first";
    CHECK(found);

    for (auto w : windows)
        delete w;
    vis::close_all();
}
