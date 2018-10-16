################################################################################
#
# Common ATDM CTest -S configuration settings for all ATDM bulids
#
# System-specific configurations are handled in the *.cmake driver files for
# each sytem in subdirectories.
#
# The CDash build name is the same as the Jenkins JOB_NAME var.
#
################################################################################

MACRO(ATDM_ASSERT_ENV_VAR_SET VAR_NAME)
  IF ("$ENV{${VAR_NAME}}" STREQUAL "")
    MESSAGE(FATAL_ERROR "Error, must set env var ${VAR_NAME}")
  ENDIF()
ENDMACRO()

ATDM_ASSERT_ENV_VAR_SET(JOB_NAME)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_KNOWN_HOSTNAME)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_KNOWN_SYSTEM_NAME)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_USE_NINJA)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_CTEST_PARALLEL_LEVEL)
ATDM_ASSERT_ENV_VAR_SET(ATDM_CONFIG_BUILD_COUNT)

SET(THIS_FILE_LIST_DIR "${CMAKE_CURRENT_LIST_DIR}")
INCLUDE("${THIS_FILE_LIST_DIR}/../../TrilinosCTestDriverCore.cmake")

INCLUDE("${THIS_FILE_LIST_DIR}/../../../std/atdm/utils/ATDMDevEnvUtils.cmake")

SET(THIS_LIST_FILE "${CMAKE_CURRENT_LIST_FILE}")

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Always assume the PWD is the root directory
  SET(CTEST_DASHBOARD_ROOT PWD)

  # Must set the Jenkins JOB_NAME which will be the CDash build name 
  SET(CTEST_BUILD_NAME "$ENV{JOB_NAME}")

  # Add this script and the shiller env script to the notes
  SET( CTEST_NOTES_FILES
    ${CTEST_NOTES_FILES}
    "${TRIBITS_PROJECT_ROOT}/cmake/std/atdm/$ENV{ATDM_CONFIG_KNOWN_SYSTEM_NAME}/environment.sh"
    "${THIS_FILE_LIST_DIR}/$ENV{ATDM_CONFIG_KNOWN_SYSTEM_NAME}/drivers/$ENV{JOB_NAME}.sh"
    "${CMAKE_CURRENT_LIST_FILE}"
    )
  
  SET(CTEST_PARALLEL_LEVEL "$ENV{ATDM_CONFIG_CTEST_PARALLEL_LEVEL}")
  
  IF ($ENV{ATDM_CONFIG_USE_NINJA})
    SET(CTEST_CMAKE_GENERATOR Ninja)
    IF ("$ENV{ATDM_CONFIG_BUILD_COUNT}" GREATER "0")
      SET(CTEST_BUILD_FLAGS "-j$ENV{ATDM_CONFIG_BUILD_COUNT} ")
    ELSE()
      SET(CTEST_BUILD_FLAGS "")
    ENDIF()
    SET(CTEST_BUILD_FLAGS "${CTEST_BUILD_FLAGS}-k 999999")
  ELSE()
    SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
    SET(CTEST_BUILD_FLAGS "-j$ENV{ATDM_CONFIG_BUILD_COUNT} -k")
  ENDIF()
  ATDM_SET_CACHE(CTEST_BUILD_FLAGS "${CTEST_BUILD_FLAGS}" CACHE STRING)
  # NOTE: Above, we need to set this as a cache var because this var is also
  # set as a cache var in ATDMDevEnvSettings.cmake that gets included below.
  
  SET(EXTRA_CONFIGURE_OPTIONS)
  
  # Must set the site name so that it does not change depending on what node
  # runs the build.
  SET(CTEST_SITE "$ENV{ATDM_CONFIG_KNOWN_HOSTNAME}")

  # See if to enable all of the packages
  MESSAGE("ENV ATDM_CONFIG_ENABLE_ALL_PACKAGE = $ENV{ATDM_CONFIG_ENABLE_ALL_PACKAGES}")
  IF ($ENV{ATDM_CONFIG_ENABLE_ALL_PACKAGES})
    SET(ATDM_ENABLE_ALL_PACKAGES TRUE)
  ELSE()
    SET(ATDM_ENABLE_ALL_PACKAGES FALSE)
  ENDIF()

  # Get set of options files (must be relative to base Trilinos dir!)
  MESSAGE("ENV ATDM_CONFIG_CONFIGURE_OPTIONS_FILES = $ENV{ATDM_CONFIG_CONFIGURE_OPTIONS_FILES}")
  IF (NOT "$ENV{ATDM_CONFIG_CONFIGURE_OPTIONS_FILES}" STREQUAL "")
    SET(ATDM_CONFIGURE_OPTIONS_FILES $ENV{ATDM_CONFIG_CONFIGURE_OPTIONS_FILES})
  ELSE()
    SET(ATDM_CONFIGURE_OPTIONS_FILES cmake/std/atdm/ATDMDevEnv.cmake)
  ENDIF()
  PRINT_VAR(ATDM_CONFIGURE_OPTIONS_FILES)

  MESSAGE("Include the configure options files at the top level to influence what package get enabled or disabled ...")
  FOREACH (CONFIG_OPTIONS_FILE ${ATDM_CONFIGURE_OPTIONS_FILES})
    SET(CONFIG_OPTIONS_FILE "${TRIBITS_PROJECT_ROOT}/${CONFIG_OPTIONS_FILE}")
    MESSAGE("Including ${CONFIG_OPTIONS_FILE} ...")
    INCLUDE("${CONFIG_OPTIONS_FILE}")
  ENDFOREACH()

  # See if "panzer" is in the job name
  STRING(REGEX MATCH "panzer" ATDM_PANZER_IN_JOB_NAME
    "${CTEST_BUILD_NAME}" )

  IF (ATDM_ENABLE_ALL_PACKAGES)
    MESSAGE("Enabling all packages by default!")
    SET(Trilinos_PACKAGES)
  ELSEIF (ATDM_PANZER_IN_JOB_NAME)
    MESSAGE("Found 'panzer' in JOB_NAME, enabling only Panzer tests")
    SET(Trilinos_PACKAGES Panzer)
  ELSE()
    # Explictly test an important subset of Trilinos packages used by the ATDM
    # APP codes that are under active development.
    SET(Trilinos_PACKAGES
      Kokkos
      Teuchos
      KokkosKernels
      Sacado
      Tpetra
      TrilinosSS
      Thyra
      Xpetra
      Zoltan2
      Belos
      Amesos2
      SEACAS
      Anasazi
      Ifpack2
      Stratimikos
      Teko
      Intrepid2
      STK
      Phalanx
      NOX
      MueLu
      Rythmos
      Tempus
      Piro
      Panzer
      )
    # NOTE: Above, we explicilty list out the packages that we want to be tested
    # which is *not* the full set of packages used by the ATDM APP codes.  The
    # other packages are not actively developed and are at less of a risk to be
    # broken.
  ENDIF()

  # Point to the ATDM Trilinos configuration
  SET(EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=${ATDM_CONFIGURE_OPTIONS_FILES}"
    "-DTrilinos_TRACE_ADD_TEST=ON"
    "-DTrilinos_ENABLE_CONFIGURE_TIMING=ON"
    )

  # Don't bother processing packages that are only implicitly enabled due to
  # enabled downstream dependencies
  SET(CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES FALSE)

  # Assume a default timeout limit of 10 min (but can be overridden in the
  # driver script and in env)
  SET_DEFAULT(CTEST_TEST_TIMEOUT 600)

  SET(Trilinos_ENABLE_CONFIGURE_TIMING ON)

  # Add this file to the list of notes files
  SET(CTEST_NOTES_FILES
    ${CTEST_NOTES_FILES}
    "${THIS_LIST_FILE}"
    )

  # Make the default branch 'develop' (but allow it to be overridden in
  # *.cmake script)
  SET_DEFAULT(Trilinos_BRANCH develop)

  # Set the default CDash Track/Group to "Specialized".  This will not trigger
  # CDash error emails for any failures.  But when the build is clean, the var
  # Trilinos_TRACK should be overridden to be "ATDM" on a build-by-build
  # basis.
  SET_DEFAULT(CTEST_TEST_TYPE Nightly)
  SET_DEFAULT(Trilinos_TRACK Specialized)

  # Don't process any extra repos
  SET(Trilinos_EXTRAREPOS_FILE NONE)

 TRILINOS_CTEST_DRIVER()

ENDMACRO()
