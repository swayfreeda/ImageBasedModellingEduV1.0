project(mvs)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")


include_directories("..")
set(HEADERS
        defines.h
        dmrecon.h
        global_view_selection.h
        image_pyramid.h
        local_view_selection.h
        mvs_tools.h
        patch_optimization.h
        patch_sampler.h
        progress.h
        settings.h
        single_view.h
        view_selection.h
        )

set(SOURCE_FILES
        dmrecon.cc
        global_view_selection.cc
        image_pyramid.cc
        local_view_selection.cc
        mvs_tools.cc
        patch_optimization.cc
        patch_sampler.cc
        single_view.cc
        )
add_library(mvs ${HEADERS} ${SOURCE_FILES})
#target_link_libraries(sfm core util features)

