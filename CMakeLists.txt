cmake_minimum_required(VERSION 3.10)

# Includes
include(CTest)

# MPI
find_package(MPI REQUIRED)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()


# Variables
## C++ Compiler Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

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

# Project Subdirectories
add_subdirectory(${SRC_DIR})
