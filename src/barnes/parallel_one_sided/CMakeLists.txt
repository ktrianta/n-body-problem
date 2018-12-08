# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE_FILES main.cpp serialization.cpp serialization.hpp)
add_executable(barnes-parallel-one-sided ${SOURCE_FILES} ${UTILS_SOURCES})
target_compile_definitions(barnes-parallel-one-sided PRIVATE MPI_IS_USED=1)
target_link_libraries(barnes-parallel-one-sided ${MPI_LIBRARIES})
set_target_properties(barnes-parallel-one-sided PROPERTIES OUTPUT_NAME prog)

install(TARGETS barnes-parallel-one-sided DESTINATION ${BARNES_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${BARNES_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${BARNES_PARALLEL_ONE_SIDED_INSTALL_BIN_DIR}/results)
