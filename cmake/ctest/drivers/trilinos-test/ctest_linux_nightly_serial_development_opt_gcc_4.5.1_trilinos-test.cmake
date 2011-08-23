
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_OPT_DEV)
SET(COMPILER_VERSION "GCC-4.5.1")
SET(ENV{LD_LIBRARY_PATH} "/home/trilinos/install/gmp-4.3.2/lib:/home/trilinos/install/mpfr2.4.2/lib:/home/trilinos/install/mpc-0.8.1/lib:/home/trilinos/gcc4.5.1/lib64:$ENV{LD_LIBRARY_PATH}")
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
#Stokhos is explicitly disabled below to prevent the package from being
#implicitly enabled.  Sundance depends on Stokhos.
SET(EXTRA_EXCLUDE_PACKAGES Phalanx Stokhos Sundance)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost_1_38_0"
  "-DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON"
  "-DCMAKE_CXX_COMPILER:FILEPATH=/home/trilinos/gcc4.5.1/bin/g++"
  "-DCMAKE_C_COMPILER:FILEPATH=/home/trilinos/gcc4.5.1/bin/gcc"
  "-DCMAKE_Fortran_COMPILER:FILEPATH=/home/trilinos/gcc4.5.1/bin/gfortran"
  "-DTrilinos_ENABLE_Stokhos:BOOL=OFF"
  "-DTPL_ENABLE_HDF5:BOOL=ON"
  "-DHDF5_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/include"
  "-DHDF5_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
