add_executable(naive-avx512 main.cpp ${UTILS_SOURCES})
target_compile_options(naive-avx512 PRIVATE -mavx -mavx2 -mavx512f -mavx512dq)
set_target_properties(naive-avx512 PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-avx512 DESTINATION ${NAIVE_AVX512_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${NAIVE_AVX512_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${NAIVE_AVX512_INSTALL_BIN_DIR}/results)
