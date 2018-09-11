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
        image.h
        image_base.h
        image_color.h
        image_drawing.h
        image_exif.h
        image_io.h
        image_tools.h
        scene.h
        view.h
        )

set(SOURCE_FILES
        bundle.cc
        bundle_io.cc
        camera.cc
        image_exif.cc
        image_io.cc
        image_tools.cc
        scene.cc
        view.cc
        )
add_library(core ${HEADERS} ${SOURCE_FILES})
target_link_libraries(core util ${PNG_LIBRARIES} ${JPEG_LIBRARIES} ${TIFF_LIBRARIES})

