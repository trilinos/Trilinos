TRIBITS_SUBPACKAGE(AdaptersSTK)

# Optional IOSS dependency
##########################
IF(${PACKAGE_NAME}_ENABLE_SEACASIoss)
   GLOBAL_SET(HAVE_IOSS ON)
ELSE()
   GLOBAL_SET(HAVE_IOSS OFF)
ENDIF()

# Optional Teko dependency
##########################
IF(${PACKAGE_NAME}_ENABLE_Teko)
   GLOBAL_SET(HAVE_TEKO ON)
ELSE()
   GLOBAL_SET(HAVE_TEKO OFF)
ENDIF()

# Optional Muelu dependency
###########################
IF(${PACKAGE_NAME}_ENABLE_MueLu)
   GLOBAL_SET(HAVE_MUELU ON)
ELSE()
   GLOBAL_SET(HAVE_MUELU OFF)
ENDIF()

# Optional Ifpack2 dependency
#############################
IF(${PACKAGE_NAME}_ENABLE_Ifpack2)
   GLOBAL_SET(HAVE_IFPACK2 ON)
ELSE()
   GLOBAL_SET(HAVE_IFPACK2 OFF)
ENDIF()

# Combined option for testing
#############################
COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_SOLVER_DRIVERS_EPETRA
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_Stratimikos ${PACKAGE_NAME}_ENABLE_Piro ${PACKAGE_NAME}_ENABLE_NOX ${PACKAGE_NAME}_ENABLE_Rythmos
  DOCSTR "Enables support for the solver drivers for the tests/examples."
    " Defaults to ON if Stratimikos, Piro, NOX, and Rythmos are"
    " enabled, otherwise defaults to OFF."
  )

# Combined option for testing
#############################
COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_STKClassic_AND_SEACAS_SUPPORT
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_STKClassic ${PACKAGE_NAME}_ENABLE_SEACASIoss ${PACKAGE_NAME}_ENABLE_SEACASExodus
  DOCSTR "Enables support for the STK mesh database and SEACAS I/O Layers."
    " Defaults to ON if STK, SEACASIoss, SEACASExodus are"
    " enabled, otherwise defaults to OFF."
  )

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()

