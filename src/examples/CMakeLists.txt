set(CALFEM_EXAMPLES
    truss3d
    frame2d
    vibration2d
    heat2d
    plane2d
    plate2d
)

# Examples using geometry and mesh generation
if(CALFEM_WITH_MESH)
    list(APPEND CALFEM_EXAMPLES mesh2d)
endif()

# Visualisation version (<example>_vis, using calfem::vis and example_vis.h) of every example
set(CALFEM_VIS_EXAMPLES)
foreach(example ${CALFEM_EXAMPLES})
    list(APPEND CALFEM_VIS_EXAMPLES ${example}_vis)
endforeach()
list(APPEND CALFEM_EXAMPLES ${CALFEM_VIS_EXAMPLES})

foreach(example ${CALFEM_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE calfem::calfem)
    set_target_properties(${example} PROPERTIES FOLDER examples)
endforeach()

# The visualisation examples show their figures in Qt windows when the Qt backend is built
if(CALFEM_WITH_QT)
    foreach(example ${CALFEM_VIS_EXAMPLES})
        target_link_libraries(${example} PRIVATE calfem::vis_qt)
        target_compile_definitions(${example} PRIVATE CALFEM_EXAMPLE_QT)
        calfem_deploy_qt_plugins(${example})
    endforeach()
endif()

# The Gmsh shared library must be next to the executables
if(CALFEM_WITH_GMSH AND Gmsh_RUNTIME)
    foreach(target ${CALFEM_EXAMPLES})
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Gmsh_RUNTIME} $<TARGET_FILE_DIR:${target}>)
    endforeach()
endif()
