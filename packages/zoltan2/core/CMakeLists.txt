IF(Zoltan2_ENABLE_ParMETIS AND NOT TPL_ENABLE_MPI)
  MESSAGE(FATAL_ERROR "ParMETIS requires that MPI is enabled. "
    "Set TPL_ENABLE_MPI ON or set TPL_ENABLE_ParMETIS OFF or set "
    "Zoltan2_ENABLE_ParMETIS OFF.")
ENDIF()

IF(Zoltan2_ENABLE_ParMETIS AND NOT HAVE_PARMETIS_VERSION_4_0_3)
  MESSAGE(FATAL_ERROR "Zoltan2 requires ParMETIS_version 4.0.3 or later. Upgrade"
    " or set TPL_ENABLE_ParMETIS OFF or set Zoltan2_ENABLE_ParMETIS OFF.")
ENDIF()

IF(Zoltan2_ENABLE_Scotch AND NOT HAVE_SCOTCH_VERSION_6_0_3)
  MESSAGE(FATAL_ERROR "Zoltan2 requires Scotch version 6.0.3 or later. Upgrade or"
    " set TPL_ENABLE_Scotch OFF or set Zoltan2_ENABLE_Scotch OFF.")
ENDIF()

# Path to test and example data directory
#KDD NOT WORKING IN NIGHTLY TEST SYSTEM
#KDD BETTER SOLUTION:  COPY THE NEEDED INPUT FILES AND USE FROM RUN DIRECTORY
#KDD SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DZ2_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test/data/\"")

# Path to Zoltan1 test directory
#KDD WILL NEED TO FIX THIS CASE AS WELL, BUT IT IS NOT USED IN NIGHTLY TESTING
#KDD SINCE NIGHTLY TESTS DO NOT ENABLE Experimental TESTS.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DZ1_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../../packages/zoltan/test/\"")

#
# A) Define the package
#

TRIBITS_PACKAGE(Zoltan2Core ENABLE_SHADOWING_WARNINGS)

#
# A2) Enable required explicit instantiation types
#     Zoltan2 requires scalar_t = int; add it here.
#     See Bug 6298.
#

FOREACH(ZNODE_TYPE  ${TpetraCore_ETI_NODES})
  FOREACH(ZLO_TYPE  ${TpetraCore_ETI_LORDS})
    FOREACH(ZGO_TYPE  ${TpetraCore_ETI_GORDS})
      TRIBITS_ADD_ETI_INSTANTIATIONS(TpetraCore
        "S={int} LO={${ZLO_TYPE}} GO={${ZGO_TYPE}} N={${ZNODE_TYPE}}")
    ENDFOREACH()
  ENDFOREACH()
ENDFOREACH()

#
# B) Set up package-specific options
#
IF (${PROJECT_NAME}_ENABLE_Epetra)
  IF(DEFINED ${PACKAGE_NAME}_ENABLE_Epetra)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_EPETRA ${${PACKAGE_NAME}_ENABLE_Epetra})
  ELSE()
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_EPETRA ON)
  ENDIF ()
ENDIF ()

IF (${PROJECT_NAME}_ENABLE_Galeri)
 GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_GALERI ON)
ENDIF ()

IF (${PROJECT_NAME}_ENABLE_Pamgen)
 GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_PAMGEN ON)
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_PURIFY
  ${PACKAGE_NAME_UC}_PURIFY
  "Enable work-arounds to purify bugs."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_UseZoltanComm
  ${PACKAGE_NAME_UC}_USEZOLTANCOMM
  "For migration in MJ, use Zoltan's comm package instead of Tpetra."
  OFF)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Machine_Dragonfly
  ${PACKAGE_NAME_UC}_MACHINE_DRAGONFLY
  "Enable Cray Dragonfly Machine representation for task mapping."
  OFF) 

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Machine_Torus
  ${PACKAGE_NAME_UC}_MACHINE_TORUS
  "Enable Torus Machine representation for task mapping."
  OFF) 

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  INCLUDE_${PACKAGE_NAME_UC}_EXPERIMENTAL
  "Enable experimental code."
  OFF)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental_Wolf
  INCLUDE_${PACKAGE_NAME_UC}_EXPERIMENTAL_WOLF
  "Enable wolf experimental code."
  OFF)

IF (${PROJECT_NAME}_ENABLE_SCOREC)
 GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_PARMA ON)
ENDIF ()

SET(HAVE_${PACKAGE_NAME_UC}_MPI ${TPL_ENABLE_MPI})

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_METIS
#  HAVE_${PACKAGE_NAME_UC}_METIS
#  "Metis is present"
#  OFF )

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PaToH
#  HAVE_${PACKAGE_NAME_UC}_PATOH
#  "PaToH is present"
#  OFF )

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ParMETIS
#  HAVE_${PACKAGE_NAME_UC}_PARMETIS
#  "ParMetis is present"
#  OFF )

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Scotch
#  HAVE_${PACKAGE_NAME_UC}_SCOTCH
#  "Scotch is present"
#  OFF )

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_OVIS
#  HAVE_${PACKAGE_NAME_UC}_OVIS
#  "Ovis is present"
#  OFF )

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_AMD
#  HAVE_${PACKAGE_NAME_UC}_AMD
#  "AMD is present"
#  OFF )

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ANASAZI
#  HAVE_${PACKAGE_NAME_UC}_ANASAZI
#  "Anasazi is present"
#  OFF )

#
# See the comment in Zoltan2_AlgScotch.hpp to learn how to
# modify the scotch library to return the maximum memory used.
#
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SCOTCH_MEMORY_REPORT
  SHOW_${PACKAGE_NAME_UC}_SCOTCH_MEMORY
  "turn on runtime display maximum memory used by scotch"
  OFF )


#INCLUDE(CheckTypeSize)
#CHECK_TYPE_SIZE("long long" HAVE_${PACKAGE_NAME_UC}_LONG_LONG)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

# This package is called Zoltan2Core from now on.
# Tests and examples belong to Zoltan2.
# TRIBITS_ADD_TEST_DIRECTORIES(test)
# TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#
#TRIBITS_EXCLUDE_FILES()


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

