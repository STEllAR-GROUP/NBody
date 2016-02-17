cmake_minimum_required(VERSION 2.8)
project(Octree CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(HPX_ROOT "$WORK/zahra/hpx")
set(HPX_DIR "${HPX_ROOT}/lib/cmake/HPX")

find_package(HPX)

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

set(SOURCE_FILES main.cpp)
add_hpx_executable(Octree
    ESSENTIAL
    SOURCES ${SOURCE_FILES}
    COMPONENT_DEPENDENCIES iostreams)
