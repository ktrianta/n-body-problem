# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(naive-parallel-avx256 main.cpp ${UTILS_SOURCES})
target_compile_options(naive-parallel-avx256 PRIVATE -mavx -mavx2)
target_compile_definitions(naive-parallel-avx256 PRIVATE MPI_IS_USED=1)
target_link_libraries(naive-parallel-avx256 ${MPI_LIBRARIES})
set_target_properties(naive-parallel-avx256 PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-parallel-avx256 DESTINATION ${NAIVE_PARALLEL_AVX256_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${NAIVE_PARALLEL_AVX256_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${NAIVE_PARALLEL_AVX256_INSTALL_BIN_DIR}/results)
