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

# Must set the Jenkins JOB_NAME which will be the CDash build name 
ATDM_ASSERT_ENV_VAR_SET(JOB_NAME)
SET(CTEST_BUILD_NAME "$ENV{JOB_NAME}")

SET(THIS_FILE_LIST_DIR "${CMAKE_CURRENT_LIST_DIR}")
INCLUDE("${THIS_FILE_LIST_DIR}/../../TrilinosCTestDriverCore.cmake")

SET(THIS_LIST_FILE "${CMAKE_CURRENT_LIST_FILE}")

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Always assume the PWD is the root directory
  SET(CTEST_DASHBOARD_ROOT PWD)

  # Point to the ATDM Trilinos configuration
  SET(EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake"
    "-DTrilinos_ENABLE_CONFIGURE_TIMING=ON"
    )

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

  # Disable the packages that are never used by the ATDM APP codes
  INCLUDE("${THIS_FILE_LIST_DIR}/../../../std/atdm/ATDMDisables.cmake")

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
