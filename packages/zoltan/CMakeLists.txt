
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Zoltan)

#
# B) Set up package-specific options
#

MESSAGE ("KDDKDD " ${CMAKE_C_COMPILER_ID} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_Fortran_COMPILER_ID})

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_OCT
  HAVE_ZOLTAN_OCT
  "Enable support for OCT in ${PACKAGE_NAME}."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_F90INTERFACE
  BUILD_ZOLTAN_F90_INTERFACE
  "Enable build of F90 interface to ${PACKAGE_NAME}."
  OFF )

OPTION(${PACKAGE_NAME}_ENABLE_CPPDRIVER
  "Enable C++ driver for ${PACKAGE_NAME}."
  ${Trilinos_ENABLE_CXX} )

IF(Trilinos_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)

  # Special Fortran compiler-specific flags
  IF (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPGI")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPGI")
    MESSAGE("KDDKDD FLAGS " ${CMAKE_C_FLAGS} " " ${CMAKE_Fortran_FLAGS})
  ENDIF()

  IF (${CMAKE_Fortran_COMPILER_ID} MATCHES "LAHEY")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DFUJITSU")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DFUJITSU")
    MESSAGE("KDDKDD FLAGS " ${CMAKE_C_FLAGS} " " ${CMAKE_Fortran_FLAGS})
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

PACKAGE_ADD_TEST_DIRECTORIES(src/driver)

IF(ZOLTAN_BUILD_ZFDRIVE)
  # CMAKE version > 2.6.4 is needed; earlier versions use C++ compiler to link
  # fdriver, while fdriver needs F90.
  PACKAGE_ADD_TEST_DIRECTORIES(src/fdriver)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
