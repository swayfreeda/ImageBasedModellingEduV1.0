project(core)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")


# find linpng
find_package(PNG REQUIRED)
if(PNG_FOUND)
    message("PNG found: ${PNG_LIBRARIES}")
    include_directories(${PNG_INCLUDE_DIRS})
    add_definitions(${PNG_DEFINITIONS})
endif()

# find libjpeg
find_package(JPEG REQUIRED)
if(JPEG_FOUND)
    message("JPEG found: ${JPEG_LIBRARIES}")
    include_directories(${JPEG_INCLUDE_DIR})
endif()

# find libtiff
find_package(TIFF REQUIRED)
if(TIFF_FOUND)
    message("TIFF found: ${TIFF_LIBRARIES}")
    include_directories(${TIFF_INCLUDE_DIR})
endif()

include_directories("..")
set(HEADERS
       defines.h
        bundle.h
        bundle_io.h
        camera.h
        depthmap.h
        image.h
        image_base.h
        image_color.h
        image_drawing.h
        image_exif.h
        image_io.h
        image_tools.h
        scene.h
        view.h
        mesh.h
        mesh_info.h
        mesh_tools.h
        mesh_io_ply.h
        mesh_io_obj.h
        mesh_io_off.h
        mesh_io_npts.h
        mesh_io_smf.h
        mesh_io_pbrt.h
        )

set(SOURCE_FILES
        bundle.cc
        bundle_io.cc
        camera.cc
        depthmap.cc
        image_exif.cc
        image_io.cc
        image_tools.cc
        scene.cc
        view.cc
        mesh.cc
        mesh_info.cc
        mesh_tools.cc
        mesh_io.cc
        mesh_io_ply.cc
        mesh_io_obj.cc
        mesh_io_off.cc
        mesh_io_npts.cc
        mesh_io_smf.cc
        mesh_io_pbrt.cc
        )
add_library(core ${HEADERS} ${SOURCE_FILES})
target_link_libraries(core util ${PNG_LIBRARIES} ${JPEG_LIBRARIES} ${TIFF_LIBRARIES})

