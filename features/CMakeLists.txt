project(features)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

#include_directories(../math)
#include_directories(../util)
include_directories("..")
set(HEADERS
        sift.h
        surf.h
        nearest_neighbor.h
        matching_base.h
        matching.h
        exhaustive_matching.h
        cascade_hashing.h
        )

set(SOURCE_FILES
        sift.cc
        surf.cc
        nearest_neighbor.cc
        matching.cc
        exhaustive_matching.cc
        cascade_hashing.cc

        )
add_library(${PROJECT_NAME} ${HEADERS} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} util core)

