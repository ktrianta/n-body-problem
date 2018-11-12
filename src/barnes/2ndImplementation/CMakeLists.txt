# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE_FILES main.cpp serialization.cpp serialization.hpp)
add_executable(barnes-2ndImplementation ${SOURCE_FILES} ${UTILS_SOURCES})
target_link_libraries(barnes-2ndImplementation ${MPI_LIBRARIES})
set_target_properties(barnes-2ndImplementation PROPERTIES OUTPUT_NAME prog)

install(TARGETS barnes-2ndImplementation DESTINATION ${BARNES_2NDIMPLEMENTATION_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${BARNES_2NDIMPLEMENTATION_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${BARNES_2NDIMPLEMENTATION_INSTALL_BIN_DIR}/results)
