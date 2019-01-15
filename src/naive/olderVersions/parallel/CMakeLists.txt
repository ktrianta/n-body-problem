# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(naive-parallel main.cpp ${UTILS_SOURCES})
target_compile_definitions(naive-parallel PRIVATE MPI_IS_USED=1)
target_link_libraries(naive-parallel ${MPI_LIBRARIES})
set_target_properties(naive-parallel PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-parallel DESTINATION ${NAIVE_PARALLEL_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${NAIVE_PARALLEL_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${NAIVE_PARALLEL_INSTALL_BIN_DIR}/results)
