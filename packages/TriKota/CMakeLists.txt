INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(TriKota DISABLE_STRONG_WARNINGS)

#
# B) Set up package-specific options
#

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
#  TRIKOTA_DEBUG
#  "Build with Debug mode."
#  OFF )

SET(TriKota_ENABLE_DakotaCMake OFF CACHE BOOL 
"Compile Dakota using CMake instead of the default autotools")

#
# C) Add the libraries, tests, and examples
#

# New support for building Dakota with CMake 
# This works by making Dakota appear as a subdirectory of TriKota
# and is enabled by configuring with TriKota_ENABLE_DakotaCMake on.
# This gives you complete control of the Dakota configure through your
# Trilinos configure scripts, makes it much easier to rebuild Dakota
# after updates, and is much faster.
#
# Currently Dakota does its own checking for libraries such as blas,
# lapack, boost, mpi, etc... which ultimately we would want to pass in.
# Also the following packages likely need to be disabled via 
# -D HAVE_X=OFF, where X is:
#   ACRO, HOPSPACK, AMPL
# Acro isn't cmake'd yet anyway, but boost_signals will be a problem.
# HOPSPACK doesn't link correctly in static builds.
# Turning AMPL on caused very strange runtime errors.
# I don't know how to turn these off by default.
IF(TriKota_ENABLE_DakotaCMake)

  # Check for Dakota source
  IF (NOT EXISTS ${PACKAGE_SOURCE_DIR}/Dakota/CMakeLists.txt)
    MESSAGE(FATAL_ERROR "\nTriKota Fatal Error: Dakota needs to 
    be untarred into the Trikota Directory before running cmake!")
  ENDIF()
  
  # DAKOTA should follow Trilinos MPI configuration
  IF(TPL_ENABLE_MPI)
    SET(DAKOTA_HAVE_MPI ON CACHE BOOL "DAKOTA MPI set by Trilinos")
  ENDIF()

  # Tell DAKOTA (and all sub-packages) to use PROJECT_NAME, e.g.,
  # Trilinos or VERA for export targets.  For a non TriBITS system
  # would instead want CMAKE_PROJECT_NAME
  SET(ExportTarget ${PROJECT_NAME} CACHE STRING 
    "Name for the export target for ${PROJECT_NAME}"
    )
  message(STATUS "TriKota setting ExportTarget=${ExportTarget}")    

  # Tell Dakota it is building inside Trilinos
  # (affects use of Teuchos and exported lists of includes and libs)
  SET(BUILD_IN_TRILINOS ON)

  # Work around what appears to be a limitation in FindBoost.cmake and
  # minimize chance Trilinos and DAKOTA pick up different Boost versions
  # (Boost_INCLUDE_DIRS doesn't seem to work as expected)
  # This may fail is Boost_INCLUDE_DIRS contains a list of DIRS!
  # Could consider mapping BoostLib_LIBRARY_DIRS to Boost_LIBRARY_DIRS as well
  IF(TPL_ENABLE_Boost)
    IF(Boost_INCLUDE_DIRS)
      SET(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIRS} CACHE FILEPATH 
	"Boost include directory set by TriKota")
    ENDIF()
  ENDIF()

  # Turn off a bunch of Dakota non-TriBITS-complient tests when buliding
  # Dakota under Trilinos/TriKota.  Until these tests can be prefixed with
  # TriKota_ and given the correct label there is no reason to enale them
  # because they will not be run in automated testing.  However, allow users
  # to enable them if they would really like.
  SET(DAKOTA_ENABLE_TESTS OFF CACHE BOOL "")
  SET(DAKOTA_ENABLE_TPL_TESTS OFF CACHE BOOL "")
  SET(HAVE_ACRO OFF CACHE BOOL "")
  SET(HAVE_AMPL OFF CACHE BOOL "")
  SET(HAVE_X_GRAPHICS OFF CACHE BOOL "")
  SET(HAVE_HOPSPACK OFF CACHE BOOL "")

  # Instruct DAKOTA to install its non-standard content to alternate locations
  SET(DAKOTA_EXAMPLES_INSTALL example/dakota CACHE FILEPATH
      "Installation destination for DAKOTA examples/ dir")
  SET(DAKOTA_TEST_INSTALL example/dakota CACHE FILEPATH
      "Installation destination for DAKOTA test/ dir")
  SET(DAKOTA_TOPFILES_INSTALL share/dakota/stable CACHE FILEPATH
      "Installation destination for DAKOTA top-level files")
  
  # Configure Dakota and its packages
  MESSAGE("-->TriKota: Configuring Dakota using CMake")
  ADD_SUBDIRECTORY(Dakota)
  MESSAGE("<--TriKota: Finished Dakota configure.")

  # Put all of Dakota's compile line defines in TriKota_config.h which
  # expands @DAKOTA_DEFINES@
  GET_PROPERTY(DAKOTA_DEFS DIRECTORY Dakota/src PROPERTY COMPILE_DEFINITIONS)
  SET(DAKOTA_DEFINES  "")
  FOREACH(DAKOTA_DEF IN LISTS DAKOTA_DEFS)
    SET(DAKOTA_DEFINES "${DAKOTA_DEFINES}#define ${DAKOTA_DEF}\n")
  ENDFOREACH()
ENDIF()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
