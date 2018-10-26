cmake_minimum_required(VERSION 3.10)

project(nbody)


set(UTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils)
set(UTILS_SOURCES ${UTILS_DIR}/initialization.cpp)
set(UTILS_HEADERS ${UTILS_DIR}/initialiaztion.hpp ${UTILS_DIR}/types.hpp)

include_directories(${UTILS_DIR})

add_subdirectory(naive/sequential)
