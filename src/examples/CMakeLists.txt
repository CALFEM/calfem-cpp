set(CALFEM_EXAMPLES
    truss3d
    frame2d
    frame2d_vis
    vibration2d
    heat2d
    plane2d
    plate2d
)

# Examples using geometry and mesh generation
if(CALFEM_WITH_MESH)
    list(APPEND CALFEM_EXAMPLES mesh2d)
endif()

foreach(example ${CALFEM_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE calfem::calfem)
    set_target_properties(${example} PROPERTIES FOLDER examples)
endforeach()

# frame2d_vis shows its figures in Qt windows when the Qt backend is built
if(CALFEM_WITH_QT)
    target_link_libraries(frame2d_vis PRIVATE calfem::vis_qt)
    target_compile_definitions(frame2d_vis PRIVATE CALFEM_EXAMPLE_QT)
endif()

# The Gmsh shared library must be next to the executables
if(CALFEM_WITH_GMSH AND Gmsh_RUNTIME)
    foreach(target ${CALFEM_EXAMPLES})
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Gmsh_RUNTIME} $<TARGET_FILE_DIR:${target}>)
    endforeach()
endif()
