project(nbody)


set(UTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils)
set(UTILS_SOURCES ${UTILS_DIR}/initialization.cpp
                  ${UTILS_DIR}/io.cpp
                  ${UTILS_DIR}/args.cpp
                  ${UTILS_DIR}/boxComputation.cpp
                  ${UTILS_DIR}/energy.cpp
                  ${UTILS_DIR}/sort.cpp)
set(UTILS_HEADERS ${UTILS_DIR}/initialiaztion.hpp
                  ${UTILS_DIR}/types.hpp
                  ${UTILS_DIR}/io.hpp
                  ${UTILS_DIR}/args.hpp
                  ${UTILS_DIR}/boxComputation.hpp
                  ${UTILS_DIR}/energy.hpp
                  ${UTILS_DIR}/sort.hpp)

include_directories(${UTILS_DIR})

add_subdirectory(naive/sequential)
add_subdirectory(naive/avx256)
#add_subdirectory(naive/avx512)
add_subdirectory(naive/parallel)
add_subdirectory(naive/parallel_avx256)
add_subdirectory(naive/parallel_one_sided)
add_subdirectory(barnes/sequential)
add_subdirectory(barnes/parallel)
add_subdirectory(barnes/parallel_one_sided)
add_subdirectory(barnes/parallel_balanced)
add_subdirectory(barnes/2ndImplementation)
add_subdirectory(barnes/parallel_serialization)
