rm -rf CMakeCache.txt
cmake \
-D CMAKE_INSTALL_PREFIX="/home/rppawlo/install_cthulhu" \
-D Trilinos_EXTRA_REPOSITORIES="Cthulhu" \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_Cthulhu:BOOL=ON \
-D Cthulhu_ENABLE_TESTS:BOOL=ON \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
-D CMAKE_SKIP_RULE_DEPENDENCY=ON \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D CMAKE_CXX_COMPILER:FILEPATH="clang++" \
-D CMAKE_C_COMPILER:FILEPATH="clang" \
-D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/local/bin/mpif77" \
-D CMAKE_LINKER:FILEPATH="clang-ar" \
-DTrilinos_EXTRA_LINK_FLAGS="-lgfortran -lm" \
../Trilinos 
