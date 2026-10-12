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

# calfem_deploy_qt_plugins(<target>)
# Copies the Qt platform (and style) plugins next to an executable after it is built, so it finds them when run
# from bin/<config> (Qt looks in <exe dir>/platforms). Only needed on Windows; elsewhere Qt finds its plugins
# through its installation prefix.
function(calfem_deploy_qt_plugins target)
    if(NOT WIN32)
        return()
    endif()
    set(plugins
        platforms Qt6::QWindowsIntegrationPlugin
        platforms Qt6::QMinimalIntegrationPlugin
        styles Qt6::QModernWindowsStylePlugin)
    while(plugins)
        list(POP_FRONT plugins plugin_dir plugin_target)
        if(TARGET ${plugin_target})
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:${target}>/${plugin_dir}"
                COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:${plugin_target}>"
                        "$<TARGET_FILE_DIR:${target}>/${plugin_dir}"
                VERBATIM)
        endif()
    endwhile()
endfunction()
