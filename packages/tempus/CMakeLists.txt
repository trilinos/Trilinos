INCLUDE(CombinedOption)

TRIBITS_PACKAGE(Tempus)

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_impl.hpp")

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()


COMBINED_OPTION( TEMPUS_ENABLE_EPETRA_STACK
  DEP_OPTIONS_NAMES Tempus_ENABLE_Epetra
                    Tempus_ENABLE_ThyraEpetraAdapters
                    Tempus_ENABLE_Belos
                    Tempus_ENABLE_Stratimikos
                    Stratimikos_ENABLE_Belos
  DOCSTR "Enables support for Epetra."
    " Defaults to ON if the Epetra, Stratimikos, and Belos are all enabled,"
    " otherwise defaults to OFF."
  )
MESSAGE(STATUS "Tempus: Epetra Stack: ${TEMPUS_ENABLE_EPETRA_STACK}")

COMBINED_OPTION( TEMPUS_ENABLE_TPETRA_STACK
  DEP_OPTIONS_NAMES Tempus_ENABLE_Tpetra
                    Tempus_ENABLE_ThyraTpetraAdapters
                    Tempus_ENABLE_Belos
                    Tempus_ENABLE_Stratimikos
                    Stratimikos_ENABLE_Belos
  DOCSTR "Enables support for Tpetra."
    " Defaults to ON if the Tpetra, Stratimikos, and Belos are all enabled,"
    " otherwise defaults to OFF."
  )
MESSAGE(STATUS "Tempus: Tpetra Stack: ${TEMPUS_ENABLE_TPETRA_STACK}")

# For use by Tempus unit tests
GLOBAL_SET(TEMPUS_UNIT_TEST_MAIN
  "${PACKAGE_SOURCE_DIR}/test/TestUtils/Tempus_UnitTestMain.cpp")

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

TRIBITS_PACKAGE_POSTPROCESS()
