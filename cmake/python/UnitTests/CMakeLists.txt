

INCLUDE(PackageAddAdvancedTest)


IF (WIN23 AND NOT CYGWIN)
  SET(DIFF_COMMAND_STUFF "${CMAKE_COMMAND}" ARGS -E compare_files)
ELSE()
  SET(DIFF_COMMAND_STUFF diff ARGS)
ENDIF()


PACKAGE_ADD_ADVANCED_TEST( XmlHtmlCDashOutputTest
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      -D Trilinos_PACKAGES_FILE:STRING=DependencyUnitTests/TrilinosPackagesForUnitTests
      -D Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON
      -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON
      -D Trilinos_DEPS_XML_OUTPUT_FILE:STRING=TrilinosPackageDependencies.xml
      -D Trilinos_DEPS_HTML_OUTPUT_FILE:STRING=TrilinosPackageDependenciesTable.html
      -D Trilinos_CDASH_DEPS_XML_OUTPUT_FILE:STRING=CDashSubprojectDependencies.xml
      -D Trilinos_ENABLE_C:BOOL=OFF
      -D Trilinos_ENABLE_CXX:BOOL=OFF
      -D Trilinos_ENABLE_Fortran:BOOL=OFF
      ${TRILINOS_HOME_DIR}
    PASS_REGULAR_EXPRESSION "Shortcircuiting after dependency tracking ..."
  TEST_1 CMND ${DIFF_COMMAND_STUFF}
    ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependencies.gold.xml
    TrilinosPackageDependencies.xml
  TEST_2 CMND ${DIFF_COMMAND_STUFF}
    ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependenciesTable.gold.html
     TrilinosPackageDependenciesTable.html
  TEST_3 CMND ${DIFF_COMMAND_STUFF}
    ${CMAKE_CURRENT_SOURCE_DIR}/CDashSubprojectDependencies.gold.xml
    CDashSubprojectDependencies.xml
  )
