# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(naive-parallel-inplace main.cpp ${UTILS_SOURCES})
target_link_libraries(naive-parallel-inplace ${MPI_LIBRARIES})
set_target_properties(naive-parallel-inplace PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-parallel-inplace DESTINATION ${NAIVE_PARALLEL_INPLACE_INSTALL_BIN_DIR})
