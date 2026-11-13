

TRIBITS_PACKAGE(Stratimikos)

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Belos_Timers
  HAVE_STRATIMIKOS_BELOS_TIMERS
  "Add timers to Belos::MultiVecTraits for Thyra MultiVector operations."
  ON )


TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Add Belos solvers directory to include path for monolithic builds
# This allows Stratimikos adapters to include Belos headers without the solvers/ prefix
IF(${PACKAGE_NAME}_ENABLE_Belos AND NOT TPL_ENABLE_Belos)
  # Monolithic build: add the Belos source directory
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../belos/tpetra/src/solvers)
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

TRIBITS_PACKAGE_POSTPROCESS()

#
# ToDos:
#
# (*) 2008/12/10: rabartl: Add general support for generating sample
# output used by doxygen documentation.
#
