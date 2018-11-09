cmake_minimum_required(VERSION 3.10)

project(nbody)


set(UTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils)
set(UTILS_SOURCES ${UTILS_DIR}/initialization.cpp
                  ${UTILS_DIR}/io.cpp
                  ${UTILS_DIR}/args.cpp)
set(UTILS_HEADERS ${UTILS_DIR}/initialiaztion.hpp
                  ${UTILS_DIR}/types.hpp
                  ${UTILS_DIR}/io.hpp
                  ${UTILS_DIR}/args.hpp)

include_directories(${UTILS_DIR})

add_subdirectory(naive/sequential)
add_subdirectory(naive/parallel)
add_subdirectory(naive/parallel_one_sided)
add_subdirectory(barnes/sequential)
add_subdirectory(barnes/parallel)
add_subdirectory(barnes/parallel_one_sided)
add_subdirectory(barnes/2ndImplementation)
