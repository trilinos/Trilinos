

TRIBITS_PACKAGE(Stratimikos)

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Belos_Timers
  HAVE_STRATIMIKOS_BELOS_TIMERS
  "Add timers to Belos::MultiVecTraits for Thyra MultiVector operations."
  ON )


TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

ADD_SUBDIRECTORY(adapters)

ADD_SUBDIRECTORY(src)

# Documentation not available in the public repository.
IF (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/CMakeLists.txt)
  ADD_SUBDIRECTORY(doc)
ENDIF()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_PACKAGE_POSTPROCESS()

#
# ToDos:
#
# (*) 2008/12/10: rabartl: Add general support for generating sample
# output used by doxygen documentation.
#
