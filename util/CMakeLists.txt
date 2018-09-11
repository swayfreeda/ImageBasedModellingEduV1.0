project(util)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

include_directories("..")
set(HEADERS
      aligned_allocator.h
        aligned_memory.h
        arguments.h
        defines.h
        exception.h
        file_system.h
        frame_timer.h
        ini_parser.h
        logging.h
        strings.h
        system.h
        timer.h
        tokenizer.h
        )

set(SOURCE_FILES
        arguments.cc
        file_system.cc
        ini_parser.cc
        system.cc

        )
add_library(util ${HEADERS} ${SOURCE_FILES})

