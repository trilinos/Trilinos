INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Phalanx)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_EXPLICIT_TEMPLATE_INSTANTIATION
  PHX_ETI
  "Build with explicit template instation enabled."
  OFF )

# If not explicitly defined by the user, pick a default node type that
# is enabled in kokkos.  This is complicated by ETI node type
# instantiation in Tpetra.  We need to use a node instantiated by
# Tpetra for panzer to build correctly.  Below we copy the logic from
# Tpetra for choosing the default node type to be in line with tpetra.
# This is somewhat fragile until we get new ETI system in place.
IF(NOT DEFINED ${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE)
  IF (Kokkos_ENABLE_Cuda)
    SET(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE "CUDA")
  ELSE()
    IF (Kokkos_ENABLE_OpenMP)
      SET(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE "OPENMP")
    ELSE()
      IF (Kokkos_ENABLE_Serial)
        SET(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE "SERIAL")
      ELSE()
        IF (Kokkos_ENABLE_Pthread)
          SET(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE "THREAD")
        ELSE()
	  MESSAGE(FATAL_ERROR "No Kokkos execution space is enabled.")
        ENDIF()
      ENDIF()
    ENDIF()
  ENDIF()
ENDIF()

IF(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE STREQUAL "OPENMP")
  SET(PHX_KOKKOS_DEVICE_TYPE_OPENMP ON)
  MESSAGE(STATUS "Execution Space: OPENMP")
ELSEIF(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE STREQUAL "CUDA")
  SET(PHX_KOKKOS_DEVICE_TYPE_CUDA ON)
  MESSAGE(STATUS "Execution Space: CUDA")
ELSEIF(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE STREQUAL "THREAD")
  SET(PHX_KOKKOS_DEVICE_TYPE_THREAD ON)
  MESSAGE(STATUS "Execution Space: THREAD")
ELSEIF(${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE STREQUAL "SERIAL")
  SET(PHX_KOKKOS_DEVICE_TYPE_SERIAL ON)
  MESSAGE(STATUS "Execution Space: SERIAL")
ELSE()
  MESSAGE(FATAL_ERROR, " The specified device type ${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE = ${${PACKAGE_NAME}_KOKKOS_DEVICE_TYPE} is not valid!  Please set this to one of the valid device types: SERIAL,OPENMP,CUDA,THREAD.")
ENDIF()

# Index size type
IF(NOT DEFINED ${PACKAGE_NAME}_INDEX_SIZE_TYPE)
  SET(${PACKAGE_NAME}_INDEX_SIZE_TYPE "UINT")
ENDIF()

IF(${PACKAGE_NAME}_INDEX_SIZE_TYPE STREQUAL "KOKKOS")
  SET(PHX_INDEX_SIZE_TYPE_KOKKOS ON)
  MESSAGE(STATUS "Index Size Type: KOKKOS")
ELSEIF(${PACKAGE_NAME}_INDEX_SIZE_TYPE STREQUAL "INT")
  SET(PHX_INDEX_SIZE_TYPE_INT ON)
  MESSAGE(STATUS "Index Size Type: INT")
ELSEIF(${PACKAGE_NAME}_INDEX_SIZE_TYPE STREQUAL "UINT")
  SET(PHX_INDEX_SIZE_TYPE_UINT ON)
  MESSAGE(STATUS "Index Size Type: UINT")
ELSEIF(${PACKAGE_NAME}_INDEX_SIZE_TYPE STREQUAL "LONGINT")
  SET(PHX_INDEX_SIZE_TYPE_LONGINT ON)
  MESSAGE(STATUS "Index Size Type: LONGINT")
ELSEIF(${PACKAGE_NAME}_INDEX_SIZE_TYPE STREQUAL "ULONGINT")
  SET(PHX_INDEX_SIZE_TYPE_ULONGINT ON)
  MESSAGE(STATUS "Index Size Type: ULONGINT")
ELSE()
  MESSAGE(FATAL_ERROR, " The specified index size type ${PACKAGE_NAME}_INDEX_SIZE_TYPE = ${${PACKAGE_NAME}_INDEX_SIZE_TYPE} is not valid!  Please set this to one of the valid types: KOKKOS,INT,UINT,LONGINT,ULONGINT.")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PHX_DEBUG
  "Build with debug code enabled. Typically, this enables array bounds checking and runtime checking of rank accessors in multidimensional arrays."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  PHX_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors wrapped around each evaluator."
  "${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}" )

IF(HAVE_INTREPID_KOKKOSCORE AND ${PACKAGE_NAME}_ENABLE_Intrepid)
  MESSAGE(STATUS "Build MDField support for Intrepid: ON")
ELSE()
  MESSAGE(STATUS "Build MDField support for Intrepid: OFF")
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

INSTALL(PROGRAMS scripts/phalanx_create_evaluator.py 
	DESTINATION bin) 

TRIBITS_PACKAGE_POSTPROCESS()
