cmake_minimum_required(VERSION 3.5)

## MPI
find_package(MPI REQUIRED)

## C++ Compiler Flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftree-vectorize -O3 -g3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -march=native -mavx -mavx2 -ftree-vectorize -g3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -mavx -mavx2 -g3 -ftree-vectorize -fopt-info-vec-optimized")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -mavx -mavx2 -g3")

## Source Directory
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)

## Tests and Resources
set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)
set(RESOURCES_DIR ${TEST_DIR}/resources)

## Binary Directories
set(NAIVE_SEQUENTIAL_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/naive/sequential)
set(NAIVE_PARALLEL_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/naive/parallel)
set(NAIVE_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/naive/parallel_one_sided)
set(BARNES_SEQUENTIAL_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/barnes/sequential)
set(BARNES_PARALLEL_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/barnes/parallel)
set(BARNES_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/barnes/parallel_one_sided)
set(BARNES_2NDIMPLEMENTATION_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/barnes/2ndImplementation)
set(BARNES_PARALLEL_SERIALIZATION_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin/barnes/parallel_serialization)

# Project Subdirectories
add_subdirectory(${SRC_DIR})
