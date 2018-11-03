project(class7)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

include_directories("../..")
include_directories(../../3rdParty/eigen)

#set(TBB_INCLUDE_DIRS "/usr/local/Cellar/tbb/2019_U1/include" CACHE PATH "Path to where the TBB include files (tbb/) reside")
#set(TBB_LIBRARY "/usr/local/Cellar/tbb/2019_U1/lib" CACHE PATH "Path to where the TBB library files (libtbb.so etc.) reside")
#list(APPEND CMAKE_MODULE_PATH "../../cmake")
###find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)


####find_package(OpenGL REQUIRED)

set(MESH_CLEAN_SOURCES
        class7_meshclean.cc)

add_executable(task7_mesh_clean ${MESH_CLEAN_SOURCES})
target_link_libraries(task7_mesh_clean util core surface)

set(TEXTURING_SOURCES
        arguments.h
        arguments.cpp
        class7_texrecon.cpp)
add_executable(task7_texturing ${TEXTURING_SOURCES})
target_link_libraries(task7_texturing mvs util core texturing coldet mrf gco)
