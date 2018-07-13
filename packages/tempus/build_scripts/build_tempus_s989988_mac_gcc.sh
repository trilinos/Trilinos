#!/bin/bash
rm -rf CMakeCache.txt CMakeFiles

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=/Users/ccober/myOffice/ATDM/build \
    -D Trilinos_ENABLE_DEBUG=OFF \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
    -D Trilinos_ENABLE_TESTS:BOOL=OFF \
    -D Trilinos_ENABLE_Tempus:BOOL=ON \
    -D Tempus_ENABLE_TESTS:BOOL=ON \
    -D Tempus_ENABLE_EXAMPLES:BOOL=ON \
    -D Tempus_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D MPI_BASE_DIR:FILEPATH=/Users/ccober/local/openmpi \
    -D CMAKE_CXX_COMPILER:FILEPATH="mpicxx" \
    -D CMAKE_C_COMPILER:FILEPATH="mpicc" \
    -D CMAKE_CXX_FLAGS:STRING="-g -Wall" \
    -D CMAKE_C_FLAGS:STRING="-g" \
    -D Trilinos_ENABLE_Fortran:BOOL=OFF \
    -D Trilinos_ENABLE_CXX11:BOOL=ON \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
    -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
    -D CMAKE_SKIP_RULE_DEPENDENCY=ON \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=OFF \
    -D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
    -D Trilinos_DEPS_XML_OUTPUT_FILE:FILEPATH="" \
    -D Trilinos_TPL_SYSTEM_INCLUDE_DIRS:BOOL=ON \
    -D BUILD_SHARED_LIBS:BOOL=FALSE \
/Users/ccober/myOffice/Trilinos

# Flags to set for building in debug mode
#    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
#    -D CMAKE_CXX_FLAGS:STRING="-g -Wall" \
