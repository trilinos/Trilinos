
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(Zoltan)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_UINT_IDS
  UNSIGNED_INT_GLOBAL_IDS
  "typedef unsigned int ZOLTAN_ID_TYPE"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_ULONG_IDS
  UNSIGNED_LONG_GLOBAL_IDS
  "typedef unsigned long ZOLTAN_ID_TYPE"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_ULLONG_IDS
  UNSIGNED_LONG_LONG_GLOBAL_IDS
  "typedef unsigned long long ZOLTAN_ID_TYPE"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_HUND
  HAVE_ZOLTAN_HUND
  "Enable support for HUND in ${PACKAGE_NAME}."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KNUTH_HASH
  HAVE_ZOLTAN_KNUTH_HASH
  "Enable support for Knuth's hash function in ${PACKAGE_NAME}."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_F90INTERFACE
  BUILD_ZOLTAN_F90_INTERFACE
  "Enable build of F90 interface to ${PACKAGE_NAME}."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PURIFY
  HAVE_PURIFY
  "Enable work-arounds to purify bugs."
  OFF )

OPTION(${PACKAGE_NAME}_ENABLE_CPPDRIVER
  "Enable C++ driver for ${PACKAGE_NAME}."
  ${${PROJECT_NAME}_ENABLE_CXX} )

IF (${PROJECT_NAME}_ENABLE_CXX) 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK")
ENDIF()

IF(${PROJECT_NAME}_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)

  # Special Fortran compiler-specific flags
  IF (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPGI")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPGI")
  ENDIF()

  IF (${CMAKE_Fortran_COMPILER_ID} MATCHES "LAHEY")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DFUJITSU")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DFUJITSU")
  ENDIF()

  SET(ZOLTAN_BUILD_ZFDRIVE YES)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

IF (NOT TPL_ENABLE_MPI)
  ADD_SUBDIRECTORY(siMPI)
  SET(ZOLTAN_BUILD_ZFDRIVE NO)
ENDIF()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(src/driver)

IF(ZOLTAN_BUILD_ZFDRIVE)
  # CMAKE version > 2.6.4 is needed; earlier versions use C++ compiler to link
  # fdriver, while fdriver needs F90.
  TRIBITS_ADD_TEST_DIRECTORIES(src/fdriver)
ENDIF()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_EXCLUDE_FILES(
  test/ch_brack2_3
  test/ch_bug
  test/ch_degenerate
  test/ch_degenerateAA
  test/ch_drake
  test/ch_ewgt
  test/ch_grid20x19
  test/ch_hammond
  test/ch_hammond2
  test/ch_nograph
  test/ch_onedbug
  test/ch_random
  test/ch_serial
  test/ch_slac
  test/ch_vwgt
  test/ch_vwgt2
  test/hg_cage10
  test/hg_diag500_4
  test/hg_ewgt
  test/hg_felix
  test/hg_ibm03
  test/hg_ml27
  test/hg_nograph
  test/hg_vwgt
  test/nem_ti_20k
  test/nem_ti_4k
  test/misc_siefert
  test/th
  test/bin
  doc/Zoltan_html/tu_html
  src/ZoltanComponent
  src/driver_old
  src/fdriver_old
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
