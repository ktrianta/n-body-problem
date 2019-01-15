# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE_FILES main.cpp octree.cpp octree.hpp)
add_executable(barnes-parallel ${SOURCE_FILES} ${UTILS_SOURCES})
target_compile_definitions(barnes-parallel PRIVATE MPI_IS_USED=1)
target_link_libraries(barnes-parallel ${MPI_LIBRARIES})
set_target_properties(barnes-parallel PROPERTIES OUTPUT_NAME prog)

install(TARGETS barnes-parallel DESTINATION ${BARNES_PARALLEL_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${BARNES_PARALLEL_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${BARNES_PARALLEL_INSTALL_BIN_DIR}/results)
