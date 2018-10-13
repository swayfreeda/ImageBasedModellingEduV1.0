cmake_minimum_required(VERSION 3.5)
project(ImageBasedModellingEdu)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

add_subdirectory(core)
#add_subdirectory(math)
add_subdirectory(util)
add_subdirectory(features)
add_subdirectory(sfm)
add_subdirectory(mvs)
add_subdirectory(examples)
