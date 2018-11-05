# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(naive-parallel-one-sided main.cpp ${UTILS_SOURCES})
target_link_libraries(naive-parallel-one-sided ${MPI_LIBRARIES})
set_target_properties(naive-parallel-one-sided PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-parallel-one-sided DESTINATION ${NAIVE_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${NAIVE_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${NAIVE_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR}/results)
