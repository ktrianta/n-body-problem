# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE_FILES main.cpp octree.cpp octree.hpp)
add_executable(barnes-sequential ${SOURCE_FILES} ${UTILS_SOURCES})
target_link_libraries(barnes-sequential ${MPI_LIBRARIES})
set_target_properties(barnes-sequential PROPERTIES OUTPUT_NAME prog)

install(TARGETS barnes-sequential DESTINATION ${BARNES_SEQUENTIAL_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${BARNES_SEQUENTIAL_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${BARNES_SEQUENTIAL_INSTALL_BIN_DIR}/results)
