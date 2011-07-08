
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG_DEV2)
SET(COMPILER_VERSION "GCC-4.6.0")
SET(ENV{LD_LIBRARY_PATH} "/home/trilinos/install/gmp-4.3.2/lib:/home/trilinos/install/mpfr2.4.2/lib:/home/trilinos/install/mpc-0.8.1/lib:/home/trilinos/gcc4.6.0/lib64:$ENV{LD_LIBRARY_PATH}")
#SET(CTEST_TEST_TYPE EXPERIMENTAL)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
#Disabling the following packages due to issues with them with the gc 4.6 compiler
#MOOCHO there were many test failures
#Piro there was one test failure
#Amesos there was one test failure
#Stratimikos there was one test failure
#Rythmos Disabling Stratimikos caused a lot of build errors
#stokhos there was one build error
SET(EXTRA_EXCLUDE_PACKAGES MOOCHO Piro Amesos Stratimikos Rythmos Stokhos)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_CXX_COMPILER:FILEPATH=/home/trilinos/gcc4.6.0/bin/g++"
  "-DCMAKE_C_COMPILER:FILEPATH=/home/trilinos/gcc4.6.0/bin/gcc"
  "-DCMAKE_Fortran_COMPILER:FILEPATH=/home/trilinos/gcc4.6.0/bin/gfortran"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost_1_46_1"
  "-DTrilinos_ENABLE_MOOCHO:BOOL=OFF"
  "-DTrilinos_ENABLE_Piro:BOOL=OFF"
  "-DTrilinos_ENABLE_Amesos:BOOL=OFF"
  "-DTrilinos_ENABLE_Stratimikos:BOOL=OFF"
  "-DTrilinos_ENABLE_Rythmos:BOOL=OFF"
  "-DTrilinos_ENABLE_Stokhos:BOOL=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
