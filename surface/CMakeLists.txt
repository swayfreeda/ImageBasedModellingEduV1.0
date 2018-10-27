project(surface)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")


include_directories("..")
set(HEADERS
        basis_function.h
        defines.h
        hermite.h
        iso_octree.h
        iso_surface.h
        mesh_clean.h
        octree.h
        sample.h
        sample_io.h
        triangulation.h
        voxel.h
        )

set(SOURCE_FILES
        basis_function.cc
        hermite.cc
        iso_octree.cc
        iso_surface.cc
        mesh_clean.cc
        octree.cc
        sample_io.cc
        triangulation.cc
        voxel.cc
        )
add_library(surface ${HEADERS} ${SOURCE_FILES})
#target_link_libraries(sfm core util features)