# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(naive-parallel main.cpp ${UTILS_SOURCES})
target_link_libraries(naive-parallel ${MPI_LIBRARIES})
set(CMAKE_CXX_COMPILER mpicxx)
set_target_properties(naive-parallel PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-parallel DESTINATION ${NAIVE_PARALLEL_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${NAIVE_PARALLEL_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${NAIVE_PARALLEL_INSTALL_BIN_DIR}/results)
