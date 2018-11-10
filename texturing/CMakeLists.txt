project(texturing)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

add_definitions(-DRESEARCH)

include_directories(..)
include_directories(../3rdParty/mrf)
include_directories(../3rdParty/coldet/src)

set(HEADERS
        defines.h
        debug.h
        util.h
        timer.h
        texturing.h
        histogram.h
        progress_counter.h
        material_lib.h
        obj_model.h
        poisson_blending.h
        rect.h
        rectangular_bin.h
        seam_leveling.h
        settings.h
        sparse_table.h
        texture_atlas.h
        texture_patch.h
        texture_view.h
        tri.h
        uni_graph.h
        timer.h
        )

set(SOURCE_FILES
        build_adjacency_graph.cpp
        build_obj_model.cpp
        calculate_data_costs.cpp
        generate_debug_embeddings.cpp
        generate_texture_atlases.cpp
        generate_texture_patches.cpp
        generate_texture_views.cpp
        global_seam_leveling.cpp
        histogram.cpp
        local_seam_leveling.cpp
        material_lib.cpp
        obj_model.cpp
        poisson_blending.cpp
        prepare_mesh.cpp
        prepare_mesh.cpp
        rectangular_bin.cpp
        seam_leveling.cpp
        texture_atlas.cpp
        texture_patch.cpp
        texture_view.cpp
        tri.cpp
        uni_graph.cpp
        view_selection.cpp
        timer.cpp
        )
add_library(texturing ${HEADERS} ${SOURCE_FILES})
#target_link_libraries(sfm core util features)

#file (GLOB HEADERS "*.h")
#file (GLOB SOURCES "[^_]*.cpp")
#
#set(LIBRARY tex)
#add_library(${LIBRARY} STATIC ${SOURCES})
#set_property(TARGET ${LIBRARY} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
#add_dependencies(${LIBRARY} ext_mve ext_rayint ext_eigen ext_mapmap)
#target_link_libraries(${LIBRARY} ${TBB_LIBRARIES} -lmve -lmve_util ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES})
#install(TARGETS ${LIBRARY} ARCHIVE DESTINATION lib)
