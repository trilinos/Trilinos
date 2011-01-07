
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosVersion.cmake")

#
# Set the options specific to this build case
#

SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME "ZOLTAN_INT_IDS")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})
#SET(CTEST_TEST_TIMEOUT 900)

#
# Set the rest of the system-specific options and run the dashboard build/test
# Platform/compiler specific options for octopi
#

# Base of Trilinos/cmake/ctest then BUILD_DIR_NAME

SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )

SET( CTEST_BUILD_FLAGS "-j8 -i" )

SET_DEFAULT( CTEST_PARALLEL_LEVEL "1" )
SET_DEFAULT(COMPILER_VERSION "GCC-3.4.6")


SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

# Only turn on PyTrilinos for shared libraries
SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES TriKota Optika)

# Purify info
SET(KDD_PURIFY "/usr/local/rational/releases/PurifyPlus.7.0/i386_linux2/bin/purify")
SET(KDD_PURIFY_ARGS "-best-effort -follow-child-processes=yes -cache-dir=/tmp/purify -chain-length=20  -windows=no ")

# Compilers that work with purify
SET(KDD_GCC  "/usr/bin/gcc346")
SET(KDD_GCXX "/usr/bin/g++346")

# Output of "mpicc --showme:compile" and "mpiCC --showme:compile"
SET(KDD_CFLAGS   "-std=c99 -DZOLTAN_ID_TYPE_UINT -m64 -g -I/opt/lam714-gcc346-pure/include -pthread")
SET(KDD_CXXFLAGS "-DZOLTAN_ID_TYPE_UINT -m64 -g -I/opt/lam714-gcc346-pure/include -pthread")

# Output of "mpiCC --showme:link"
SET(KDD_LINKFLAGS "-m64 -L/opt/lam714-gcc346-pure/lib -llammpio -llammpi++ -llamf77mpi -lmpi -llam -laio -laio -lutil -ldl")

# MPI info; needed for mpirun; also need this in path.
set(KDD_MPI "/opt/lam714-gcc346-pure/bin")

set(TDD_HTTP_PROXY "http://sonproxy.sandia.gov:80/")

SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DTPL_ENABLE_MPI:BOOL=ON"
    "-DMPI_USE_COMPILER_WRAPPERS:BOOL=OFF"
    "-DMPI_BIN_DIR:STRING=${KDD_MPI}"
    "-DCMAKE_C_COMPILER:STRING=${KDD_PURIFY}" 
    "-DCMAKE_C_FLAGS:STRING=${KDD_PURIFY_ARGS} ${KDD_GCC} ${KDD_CFLAGS}"
    "-DCMAKE_CXX_COMPILER:STRING=${KDD_PURIFY}" 
    "-DCMAKE_CXX_FLAGS:STRING=${KDD_PURIFY_ARGS} ${KDD_GCXX} ${KDD_CXXFLAGS}"
    "-DTrilinos_EXTRA_LINK_FLAGS:STRING=${KDD_LINKFLAGS}"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON "
    "-DMPI_EXEC_MAX_NUMPROCS:STRING=11 "
    "-DTrilinos_ENABLE_Fortran:BOOL=OFF "
    "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF "
    "-DTrilinos_ENABLE_EXAMPLES:BOOL=ON "
    "-DTrilinos_VERBOSE_CONFIGURE:BOOL=ON "
    "-DTrilinos_ENABLE_Zoltan:BOOL=ON "
    "-DZoltan_ENABLE_EXAMPLES:BOOL=OFF "
    "-DZoltan_ENABLE_TESTS:BOOL=ON "
    "-DZoltan_ENABLE_ParMETIS:BOOL=ON "
    "-DParMETIS_LIBRARY_DIRS:FILEPATH=/Net/local/proj/zoltan/arch/linux64/lib/lam/ParMETIS3" 
    "-DParMETIS_INCLUDE_DIRS:FILEPATH=/Net/local/proj/zoltan/arch/all/src/ParMETIS3" 
    "-DZoltan_ENABLE_Scotch:BOOL=OFF"
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DTrilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF"
  )

TRILINOS_CTEST_DRIVER()
