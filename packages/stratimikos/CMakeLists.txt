

TRIBITS_PACKAGE(Stratimikos)

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# Original behavior was to include ThyraEpetraAdapters as a required dependency.
# Modify this to make sure that ThyraEpetraAdapters is enabled if Epetra is, and
# that ThyraTpetraAdapters is enabled if Tpetra is
IF(${PACKAGE_NAME}_ENABLE_Epetra )
  SET(${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters ON)
ENDIF()
IF(${PACKAGE_NAME}_ENABLE_Tpetra )
  SET(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters ON)
ENDIF()

ADD_SUBDIRECTORY(adapters)

ADD_SUBDIRECTORY(src)

# Documentation not available in the public repository.
IF (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/CMakeLists.txt)
  ADD_SUBDIRECTORY(doc)
ENDIF()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_PACKAGE_POSTPROCESS()

#
# ToDos:
#
# (*) 2008/12/10: rabartl: Add general support for generating sample
# output used by doxygen documentation.
#
