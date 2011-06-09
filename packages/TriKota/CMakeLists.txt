INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(TriKota DISABLE_STRONG_WARNINGS)

#
# B) Set up package-specific options
#

#ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
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
  
  # Tell Dakota it is building inside Trilinos
  SET(BUILD_IN_TRILINOS ON)
  
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

PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
