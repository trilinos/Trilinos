# First copy the CMakeList.txt so we can build the test
configure_file(${PACKAGE_SOURCE_DIR}/install_test/CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/source/CMakeLists.txt)

# Second copy the source files needed to the build area
file(COPY ${PACKAGE_SOURCE_DIR}/perf_test/sparse/KokkosSparse_pcg.hpp DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/source)
file(COPY ${PACKAGE_SOURCE_DIR}/perf_test/sparse/KokkosSparse_block_pcg.cpp DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/source)

# Third write a configure file that can be invoked to test the library installation
configure_file(${PACKAGE_SOURCE_DIR}/install_test/run_install_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/run_install_test.sh @ONLY)

# Fourth create the build directory where the installation of the cg example will take place
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build)

add_test(NAME install_test COMMAND /bin/bash ${CMAKE_CURRENT_BINARY_DIR}/run_install_test.sh)
# KOKKOSKERNELS_ADD_TEST(NAME "install_test" COMMAND "${CMAKE_CURRENT_BINARY_DIR}/run_install_test.sh" COMPONENTS Sparse)
