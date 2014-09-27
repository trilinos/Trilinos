
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

TRIBITS_PACKAGE(Stratimikos)

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

ADD_SUBDIRECTORY(adapters)

ADD_SUBDIRECTORY(src)

# Documentation not available in the public repository.
IF (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc)
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
