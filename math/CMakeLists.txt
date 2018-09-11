project(math)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

include_directories("..")
set(HEADERS
      accum.h
        algo.h
        defines.h
        matrix.h
        matrix_tools.h
        vector.h
        )

add_library(math ${HEADERS})

