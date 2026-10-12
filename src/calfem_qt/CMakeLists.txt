# calfem::vis_qt - Qt backend for calfem::vis (interactive figure windows, PNG/PDF output).
# Qt is LGPL licensed; vcpkg builds it as shared libraries on Windows.

add_library(calfem_vis_qt
    ${PROJECT_SOURCE_DIR}/include/calfem/vis_qt.h
    figure_widget.h
    figure_widget.cpp
    painter.h
    painter.cpp
    qt_backend.cpp)
add_library(calfem::vis_qt ALIAS calfem_vis_qt)

target_link_libraries(calfem_vis_qt
    PUBLIC calfem::calfem Qt6::Widgets)

target_compile_features(calfem_vis_qt PUBLIC cxx_std_17)

set_target_properties(calfem_vis_qt PROPERTIES
    EXPORT_NAME vis_qt
    VERSION ${PROJECT_VERSION}
    DEBUG_POSTFIX d)

if(MSVC)
    target_compile_options(calfem_vis_qt PRIVATE /W4 /permissive-)
else()
    target_compile_options(calfem_vis_qt PRIVATE -Wall -Wextra)
endif()
