add_executable(naive-avx256 main.cpp ${UTILS_SOURCES})
target_compile_options(naive-avx256 PRIVATE -mavx -mavx2)
set_target_properties(naive-avx256 PROPERTIES OUTPUT_NAME prog)

install(TARGETS naive-avx256 DESTINATION ${NAIVE_AVX256_INSTALL_BIN_DIR})
install(DIRECTORY ${RESOURCES_DIR} DESTINATION ${NAIVE_AVX256_INSTALL_BIN_DIR})
install(DIRECTORY DESTINATION ${NAIVE_AVX256_INSTALL_BIN_DIR}/results)
