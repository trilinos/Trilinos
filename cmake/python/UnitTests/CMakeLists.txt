
INCLUDE(Multiline_Set)
INCLUDE(Generate_Script_File_Target)

FUNCTION(CREATE_DEPENCENCY_OUTPUT_TEST_CASE TEST_NAME TEST_PASS_COMMAND)

  MULTILINE_SET(TEST_CMND
    "#!/bin/bash\n"
    "if [ ! -e ${TEST_NAME} ]\n"
    "  then mkdir ${TEST_NAME}\n"
    "fi\n"
    "cd ${TEST_NAME}\n"
    "echo $PWD\n"
    "if [ -e CMakeCache.txt ]\n"
    "  then rm CMakeCache.txt\n"
    "fi\n"
    "${CMAKE_COMMAND}"
    "  -C ${TRILINOS_BUILD_DIR}/CMakeCache.txt"
    "  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
    "  -D Trilinos_PACKAGES_FILE:STRING=DependencyUnitTests/TrilinosPackagesForUnitTests"
    "  -D Trilinos_ENABLE_TESTS:BOOL=ON"
    "  -D Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON"
    "  -D Trilinos_DEPS_XML_OUTPUT_FILE:STRING=$PWD/TrilinosPackageDependencies.xml"
    "  ${TRILINOS_HOME_DIR}\n"
    "${TEST_PASS_COMMAND}\n"
    )
  
  SET(TEST_SCRIPT_FILE "${TEST_NAME}.sh")

  GENERATE_SCRIPT_FILE_TARGET(${TEST_SCRIPT_FILE} ${TEST_CMND}
    ${TRILINOS_HOME_DIR}/CMakeLists.txt
    ${TRILINOS_HOME_DIR}/cmake/python/TrilinosDependencies.py
    ${TRILINOS_HOME_DIR}/cmake/python/dump-package-dep-table.py
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    )
  
  ADD_TEST( ${TEST_NAME} ${TEST_SCRIPT_FILE} )

ENDFUNCTION()


CREATE_DEPENCENCY_OUTPUT_TEST_CASE(TrilinosOutputDep_XmlOutputTest
  "diff -s ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependencies.gold.xml TrilinosPackageDependencies.xml" )


CREATE_DEPENCENCY_OUTPUT_TEST_CASE(TrilinosOutputDep_HtmlOutputTest
  "diff -s ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependenciesTable.gold.html TrilinosPackageDependenciesTable.html" )
