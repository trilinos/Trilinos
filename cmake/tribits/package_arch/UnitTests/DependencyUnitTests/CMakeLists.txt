# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER


INCLUDE(TribitsAddAdvancedTest)

INCLUDE(ParseVariableArguments)

# MockTrilinos project
SET(MOCK_PROJECT_NAME Trilinos)
SET(MOCK_TRILINOS_DIR
  ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_MOCK_TRILINOS_DIR})


#####################################################################
#
# Unit tests for just creating the depencency XML files
#
#####################################################################


FUNCTION(CREATE_DEPENDENCY_XML_DUMP_TEST_CASE TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "ARGS;PASS_REGULAR_EXPRESSION_ALL;COMPARE_XML_FILE"
     #options
     ""
     ${ARGN}
     )

  SET(EXTENDED_TEST_NAME DepXmlDumpTests_${TEST_NAME})
  SET(FULL_TEST_NAME ${PACKAGE_NAME}_${EXTENDED_TEST_NAME})

  SET(TEST_OUTPUT_REL_DIR ${FULL_TEST_NAME})
  SET(TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_OUTPUT_REL_DIR})

  TRIBITS_ADD_ADVANCED_TEST( ${EXTENDED_TEST_NAME}
    OVERALL_WORKING_DIRECTORY ${TEST_OUTPUT_REL_DIR}
    OVERALL_NUM_MPI_PROCS 1
    TEST_0 CMND "${CMAKE_COMMAND}"
      ARGS
        -DPROJECT_SOURCE_DIR:PATH=${MOCK_TRILINOS_DIR}
        -D${MOCK_PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -D${MOCK_PROJECT_NAME}_NATIVE_REPOSITORIES=.
        -D${MOCK_PROJECT_NAME}_DEPS_XML_OUTPUT_FILE=${TEST_OUTPUT_DIR}/${MOCK_PROJECT_NAME}PackageDependencies.xml
        -D${MOCK_PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE=CDashSubprojectDependencies.xml
        -D${MOCK_PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE=${MOCK_PROJECT_NAME}PackageDependenciesTable.html
        ${PARSE_ARGS}
        -P ${${PROJECT_NAME}_TRIBITS_DIR}/package_arch/TribitsDumpDepsXmlScript.cmake
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    TEST_1 CMND "${CMAKE_COMMAND}"
      ARGS
        -E compare_files
        ${TEST_OUTPUT_DIR}/${MOCK_PROJECT_NAME}PackageDependencies.xml
        ${${PROJECT_NAME}_TRIBITS_DIR}/python/UnitTests/${PARSE_COMPARE_XML_FILE}
    XHOSTTYPE Windows
    )

ENDFUNCTION()


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  BasePackagesTpls
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='24'"
    "Trilinos_NUM_TPLS='19'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithPreCopyrightTrilinos
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=preCopyrightTrilinos
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithPreCopyrightTrilinosExtraTrilinosRepo
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=preCopyrightTrilinos,extraTrilinosRepo
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='22'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.preCopyrightTrilinos.extraTrilinosRepo.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithExtraRepoOnePackage
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=extraRepoOnePackage
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='24'"
    "Trilinos_NUM_TPLS='19'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.extraRepoOnePackage.gold.xml
  )


#####################################################################
#
# Unit tests for dependency handling for reduced set of packages
#
#####################################################################

# Reduced MockTrilinos project
SET(REDUCED_MOCK_PROJECT_NAME Trilinos)
SET(REDUCED_MOCK_TRILINOS_DIR
  ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_REDUCED_MOCK_TRILINOS_DIR})

# Expected output files directory
SET(REDUCED_MOCK_TRILINOS_OUTPUT_DIR
  "${${PROJECT_NAME}_TRIBITS_DIR}/doc/developers_guide/ReducedMockTrilinosOutput")

#
# NOTE: These tests compare against expected output stored in files under the
# in the developers_guide/ReducedMockTrilinosOuptut/ directory and are used in
# the TribitsDevelopersGuide.rst document.  Therefore, if the behavior of
# TriBITS changes in such a way as to break these tests, then if any of the
# contents of the files under this directory need to change, then one must
# also look at the TribitsDevelopersGuilde.rst document and update it for the
# new behvior of TriBITs!
#


FUNCTION(CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "ARGS;PASS_REGULAR_EXPRESSION_ALL"
     #options
     ""
     ${ARGN}
     )

  TRIBITS_ADD_ADVANCED_TEST( DepTestsReduced_${TEST_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    TEST_0 CMND "${CMAKE_COMMAND}"
      ARGS
        ${PARSE_ARGS}
        ${REDUCED_MOCK_TRILINOS_DIR}
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    )

ENDFUNCTION()


#
# Read in an expected output file and create a REGEX array to be used in
# checking test results.
#

FUNCTION(READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY  EXPECTED_OUTPUT_FILE
  EXPECTED_OUTPUT_REGEX_ARRAY_OUT
  )
  FILE(READ ${EXPECTED_OUTPUT_FILE} EXPECTED_RAW_STR_ALL)
  #PRINT_VAR(EXPECTED_RAW_STR_ALL)
  SPLIT("${EXPECTED_RAW_STR_ALL}" "\n" EXPECTED_RAW_STR_ARRAY)
  #PRINT_VAR(EXPECTED_RAW_STR_ARRAY)
  SET(EXPECTED_REGEX_ARRAY)
  FOREACH(EXPECTED_RAW_STR ${EXPECTED_RAW_STR_ARRAY})
    #PRINT_VAR(EXPECTED_RAW_STR)
    IF (NOT "${EXPECTED_RAW_STR}" STREQUAL "")
      SET(EXPECTED_REGEX_STR "${EXPECTED_RAW_STR}")
      STRING(REPLACE "(" "." EXPECTED_REGEX_STR "${EXPECTED_REGEX_STR}") 
      STRING(REPLACE ")" "." EXPECTED_REGEX_STR "${EXPECTED_REGEX_STR}") 
      STRING(REPLACE "[" "." EXPECTED_REGEX_STR "${EXPECTED_REGEX_STR}") 
      STRING(REPLACE "]" "." EXPECTED_REGEX_STR "${EXPECTED_REGEX_STR}") 
      STRING(REPLACE "*" "." EXPECTED_REGEX_STR "${EXPECTED_REGEX_STR}") 
      LIST(APPEND EXPECTED_REGEX_ARRAY "${EXPECTED_REGEX_STR}")
      #PRINT_VAR(EXPECTED_REGEX_ARRAY)
    ENDIF()
  ENDFOREACH()
  SET(${EXPECTED_OUTPUT_REGEX_ARRAY_OUT} ${EXPECTED_REGEX_ARRAY} PARENT_SCOPE)
ENDFUNCTION()



#
# NoEnables
#

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/NoEnables.txt"
  RMT_NoEnables_REGEX_ARRAY
  )
#PRINT_VAR(REDUCED_EXPECTED_DEPENDENICES_REGEX_ARRAY)


CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  NoEnables
  ARGS
  PASS_REGULAR_EXPRESSION_ALL
    ${RMT_NoEnables_REGEX_ARRAY}
  )


#
# EnableThyra_EnableTests
#

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/EnableThyra_EnableTests.txt"
  EnableThyra_EnableTests_REGEX_ARRAY
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_EnableTests
  ARGS
    -DTrilinos_ENABLE_Thyra:BOOL=ON
    -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    ${EnableThyra_EnableTests_REGEX_ARRAY}
  )


#
# EnableThyra_EnableTests_EnableBoost_ST
#

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/EnableThyra_EnableTests_EnableBoost_ST.txt"
  EnableThyra_EnableTests_EnableBoost_ST_REGEX_ARRAY
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_EnableTests_EnableBoost_ST
  ARGS
    -DTPL_ENABLE_Boost:BOOL=ON
    -DTrilinos_ENABLE_Thyra:BOOL=ON
    -DTrilinos_ENABLE_TESTS:BOOL=ON
    -DTrilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    ${EnableThyra_EnableTests_EnableBoost_ST_REGEX_ARRAY}
  )


#
# EnableThyra_EnableTests_EnableBoost_ST_RequiredOnly
#

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/EnableThyra_EnableTests_EnableBoost_ST_RequiredOnly.txt"
  EnableThyra_EnableTests_EnableBoost_ST_RequiredOnly_REGEX_ARRAY
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_EnableTests_EnableBoost_ST_RequiredOnly
  ARGS
    -DTPL_ENABLE_Boost:BOOL=ON
    -DTrilinos_ENABLE_Thyra:BOOL=ON
    -DTrilinos_ENABLE_TESTS:BOOL=ON
    -DTrilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    ${EnableThyra_EnableTests_EnableBoost_ST_RequiredOnly_REGEX_ARRAY}
  )


#
# EnableThyra_DisableEpetra_EnableTests
#

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/EnableThyra_DisableEpetra_EnableTests.txt"
  EnableThyra_DisableEpetra_EnableTests_REGEX_ARRAY
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_DisableEpetra_EnableTests
  ARGS
    -DTrilinos_ENABLE_Thyra:BOOL=ON
    -DTrilinos_ENABLE_Epetra:BOOL=OFF
    -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    ${EnableThyra_DisableEpetra_EnableTests_REGEX_ARRAY}
  )

#
# EnableEpetra_EnableRTOp_DisableTeuchos
#

SET(TEST_NAME "EnableEpetra_EnableRTOp_DisableTeuchos")

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/${TEST_NAME}.txt"
  ${TEST_NAME}_REGEX_ARRAY
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  ${TEST_NAME}
  ARGS
    -DTrilinos_ENABLE_Epetra:BOOL=ON
    -DTrilinos_ENABLE_RTOp:BOOL=ON
    -DTrilinos_ENABLE_Teuchos:BOOL=OFF
   PASS_REGULAR_EXPRESSION_ALL
    ${${TEST_NAME}_REGEX_ARRAY}
  )


#
# EnableEpetra_EnableRTOp_DisableTeuchos_DisableEnabledFwdDepPackages
#

SET(TEST_NAME "EnableEpetra_EnableRTOp_DisableTeuchos_DisableEnabledFwdDepPackages")

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/${TEST_NAME}.txt"
  ${TEST_NAME}_REGEX_ARRAY
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  ${TEST_NAME}
  ARGS
    -DTrilinos_ENABLE_Epetra:BOOL=ON
    -DTrilinos_ENABLE_RTOp:BOOL=ON
    -DTrilinos_ENABLE_Teuchos:BOOL=OFF
    -DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
   PASS_REGULAR_EXPRESSION_ALL
    ${${TEST_NAME}_REGEX_ARRAY}
  )


#
# EnableAllPackages_DumpDependencies
#

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/ExpectedDependencies.txt"
  REDUCED_EXPECTED_DEPENDENICES_REGEX_ARRAY
  )
#PRINT_VAR(REDUCED_EXPECTED_DEPENDENICES_REGEX_ARRAY)

READ_EXPECTED_OUTPUT_FILE_INTO_REGEX_ARRAY(
  "${REDUCED_MOCK_TRILINOS_OUTPUT_DIR}/EnableAllPackages.txt"
  REDUCED_EnableAllPackages_REGEX_ARRAY
  )
#PRINT_VAR(REDUCED_EnableAllPackages_REGEX_ARRAY)

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_DumpDependencies
  ARGS
    -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
    -DTrilinos_DUMP_PACKAGE_DEPENDENCIES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    ${REDUCED_EXPECTED_DEPENDENICES_REGEX_ARRAY}
    ${REDUCED_EnableAllPackages_REGEX_ARRAY}
    "Final set of enabled packages:  Teuchos RTOp Epetra Thyra 4"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra ThyraCoreLibs ThyraEpetra Thyra 6"
    "Final set of non-enabled packages:  Triutils EpetraExt 2"
    "Final set of non-enabled SE packages:  Triutils EpetraExt ThyraGoodStuff ThyraCrazyStuff ThyraEpetraExt 5"
  )

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_VerboseConfigure
  ARGS
    -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
    -DTrilinos_VERBOSE_CONFIGURE:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Teuchos_FULL_ENABLED_DEP_PACKAGES=''"
    "RTOp_FULL_ENABLED_DEP_PACKAGES='Teuchos'"
    "Epetra_FULL_ENABLED_DEP_PACKAGES=''"
    "ThyraCoreLibs_FULL_ENABLED_DEP_PACKAGES='RTOp.Teuchos'"
    "ThyraEpetra_FULL_ENABLED_DEP_PACKAGES='ThyraCoreLibs.Epetra.RTOp.Teuchos'"
    "Thyra_FULL_ENABLED_DEP_PACKAGES='ThyraEpetra.ThyraCoreLibs.Epetra.RTOp.Teuchos'"
  )


#
# EnableAllPackages_RequiredOnly
#

CREATE_REDUCED_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_RequiredOnly
  ARGS
    -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Thyra 4"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra ThyraCoreLibs ThyraEpetra Thyra 6"
    "Final set of non-enabled packages:  Triutils EpetraExt 2"
    "Final set of non-enabled SE packages:  Triutils EpetraExt ThyraGoodStuff ThyraCrazyStuff ThyraEpetraExt "
  )


#####################################################################
#
# Unit tests for dependency handling for full set of packages
#
#####################################################################


FUNCTION(CREATE_DEPENDENCY_HANDLING_TEST_CASE TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "ARGS;PASS_REGULAR_EXPRESSION_ALL"
     #options
     ""
     ${ARGN}
     )

  TRIBITS_ADD_ADVANCED_TEST( DepTests_${TEST_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    TEST_0 CMND "${CMAKE_COMMAND}"
      ARGS
        -D${MOCK_PROJECT_NAME}_TRIBITS_DIR:FILEPATH=${${PROJECT_NAME}_TRIBITS_DIR}
        -D${MOCK_PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON
        -D${MOCK_PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:STRING=
        -D${MOCK_PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE:STRING=
        -D${MOCK_PROJECT_NAME}_ENABLE_C:BOOL=OFF
        -D${MOCK_PROJECT_NAME}_ENABLE_CXX:BOOL=OFF
        -D${MOCK_PROJECT_NAME}_ENABLE_Fortran:BOOL=OFF
        ${PARSE_ARGS}
        ${MOCK_TRILINOS_DIR}
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    )

ENDFUNCTION()


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "WARNING: MissingPackage is being ignored since its directory is missing and MissingPackage_ALLOW_MISSING_EXTERNAL_PACKAGE = TRUE"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableTests
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "Setting Teuchos_ENABLE_TESTS=ON"
    "Setting RTOp_ENABLE_TESTS=ON"
    "Setting Epetra_ENABLE_TESTS=ON"
    "Setting ML_ENABLE_TESTS=ON"
    "Setting Stratimikos_ENABLE_TESTS=ON"
    "Setting RBGen_ENABLE_TESTS=ON"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableEpetra_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Zoltan Shards Tpetra Thyra Intrepid ML RBGen 10"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Teuchos:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_RTOp=OFF because RTOp has a required library dependence on disabled package Teuchos"
    "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableEpetra_DisableTeuchos_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Zoltan Shards ML 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Teuchos:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllOptionalPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_AztecOO:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Not Setting Trilinos_ENABLE_Belos=ON even though Stratimikos has an optional dependence on Belos because Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even though Intrepid has an optional dependence on Sacado because Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 17"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Amesos:BOOL=ON -DTrilinos_ENABLE_Zoltan:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Shards Triutils EpetraExt Thyra AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos 14"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyraCoreLibs_EnableTests
  ARGS -DTrilinos_ENABLE_ThyraCoreLibs:BOOL=ON -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_Thyra=ON because Trilinos_ENABLE_ThyraCoreLibs=ON"
    "Setting Thyra_ENABLE_TESTS=ON because ThyraCoreLibs_ENABLE_TESTS=ON"
    "Setting Thyra_ENABLE_EXAMPLES=ON because ThyraCoreLibs_ENABLE_EXAMPLES=ON"
    "Final set of enabled packages:  Teuchos RTOp Thyra 3"
    "Final set of enabled SE packages:  Teuchos RTOp ThyraCoreLibs Thyra 4"    
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Thyra:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .by user.:  Thyra 1"
    "Enabling subpackages for hard enables of parent packages due to Trilinos_ENABLE_<PARENT_PACKAGE>=ON ..."
    "Setting subpackage enable Trilinos_ENABLE_ThyraCoreLibs=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetra=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetraExt=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Setting subpackage enable Trilinos_ENABLE_ThyraTpetra=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Sweep forward enabling all forward library dependent packages because Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ..."
    "Setting Trilinos_ENABLE_Stratimikos=ON because Trilinos_ENABLE_ThyraCoreLibs=ON"
    "Final set of enabled packages:  .+ Thyra .+ Stratimikos"
    "Final set of non-enabled SE packages:  .+ ThyraGoodStuff ThyraCrazyStuff"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_DisableThyraEpetraExt
  ARGS -DTrilinos_ENABLE_Thyra:BOOL=ON -DTrilinos_ENABLE_ThyraEpetraExt:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Setting subpackage enable Trilinos_ENABLE_ThyraCoreLibs=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetra=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Setting subpackage enable Trilinos_ENABLE_ThyraTpetra=ON because parent package Trilinos_ENABLE_Thyra=ON"
    "Setting Thyra_ENABLE_ThyraCrazyStuff=OFF because Thyra has an optional library dependence on disabled package ThyraCrazyStuff"
    "Setting Trilinos_ENABLE_Stratimikos=OFF because Stratimikos has a required library dependence on disabled package ThyraEpetraExt"
    "Setting Thyra_ENABLE_ThyraEpetraExt=OFF because Thyra has an optional library dependence on disabled package ThyraEpetraExt"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils Tpetra Thyra 6"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Triutils Tpetra ThyraCoreLibs ThyraEpetra ThyraTpetra Thyra 9"
    "Final set of non-enabled packages:  TrilinosFramework Zoltan Shards EpetraExt Stokhos Sacado Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 18"
    "Final set of non-enabled SE packages:  TrilinosFramework Zoltan Shards EpetraExt Stokhos Sacado ThyraGoodStuff ThyraCrazyStuff ThyraEpetraExt Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 21"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyraEpetraExt
  ARGS -DTrilinos_ENABLE_ThyraEpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .by user.:  0"
    "Explicitly enabled SE packages on input .by user.:  ThyraEpetraExt 1"
    "Setting Trilinos_ENABLE_ThyraEpetra=ON because ThyraEpetraExt has a required dependence on ThyraEpetra"
    "Setting Trilinos_ENABLE_ThyraCoreLibs=ON because ThyraEpetra has a required dependence on ThyraCoreLibs"
    "Setting Trilinos_ENABLE_Teuchos=ON because ThyraCoreLibs has a required dependence on Teuchos"
    "Setting Trilinos_ENABLE_Thyra=ON because Trilinos_ENABLE_ThyraCoreLibs=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra 7"
    "Final set of non-enabled SE packages:  .+ ThyraGoodStuff ThyraCrazyStuff ThyraTpetra .+"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableThyra_EnableThyraEpetraExt
  ARGS -DTrilinos_ENABLE_Thyra:BOOL=OFF -DTrilinos_ENABLE_ThyraEpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetraExt=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Final set of enabled packages:  0"
    "Final set of enabled SE packages:  0"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_DisableThyraEpetra_EnableThyraEpetraExt_DisableEnabledDownstreamPackages
  ARGS -DTrilinos_ENABLE_Thyra:BOOL=ON -DTrilinos_ENABLE_ThyraEpetra:BOOL=OFF
    -DTrilinos_ENABLE_ThyraEpetraExt:BOOL=ON
    -DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .by user.:  Thyra 1"
    "Explicitly enabled SE packages on input .by user.:  ThyraEpetraExt Thyra 2"
    "Explicitly disabled SE packages on input .by user or by default.:  .+ ThyraEpetra"
    "WARNING: Setting Trilinos_ENABLE_ThyraEpetraExt=OFF which was 'ON' because ThyraEpetraExt has a required library dependence on disabled package ThyraEpetra but Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON!"
    "Setting Thyra_ENABLE_ThyraEpetra=OFF because Thyra has an optional library dependence on disabled package ThyraEpetra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils Tpetra Thyra 6"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Triutils Tpetra ThyraCoreLibs ThyraTpetra Thyra 8"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_DisableEpetra_EnableThyraEpetraExt_DisableEnabledDownstreamPackages
  ARGS -DTrilinos_ENABLE_Thyra:BOOL=ON -DTrilinos_ENABLE_Epetra:BOOL=OFF
    -DTrilinos_ENABLE_ThyraEpetraExt:BOOL=ON 
    -DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .by user.:  Thyra 1"
    "Explicitly enabled SE packages on input .by user.:  ThyraEpetraExt Thyra 2"
    "Explicitly disabled packages on input .by user or by default.:  Epetra"
    "Setting Trilinos_ENABLE_ThyraEpetra=OFF because ThyraEpetra has a required library dependence on disabled package Epetra"
    "WARNING: Setting Trilinos_ENABLE_ThyraEpetraExt=OFF which was 'ON' because ThyraEpetraExt has a required library dependence on disabled package EpetraExt but Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON!"
    "Final set of enabled packages:  Teuchos RTOp Tpetra Thyra 4"
    "Final set of enabled SE packages:  Teuchos RTOp Tpetra ThyraCoreLibs ThyraTpetra Thyra 6"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_DisableThyra
  ARGS -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Thyra=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Setting subpackage enable Trilinos_ENABLE_ThyraCoreLibs=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraGoodStuff=OFF because parent package Trilinos_ENABLE_Thyra=OFF" 
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetra=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetraExt=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraTpetra=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "ERROR: Setting Trilinos_ENABLE_Stratimikos=OFF which was 'ON' because Stratimikos has a required library dependence on disabled package ThyraCoreLibs"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_DisableThyra_DisableEnabledDownstreamPackages
  ARGS -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Thyra=OFF
    -DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting subpackage enable Trilinos_ENABLE_ThyraCoreLibs=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraGoodStuff=OFF because parent package Trilinos_ENABLE_Thyra=OFF" 
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetra=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetraExt=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraTpetra=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "WARNING: Setting Trilinos_ENABLE_Stratimikos=OFF which was 'ON' because Stratimikos has a required library dependence on disabled package ThyraCoreLibs but Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON"
    "Final set of enabled packages:  0"
    "Final set of non-enabled packages:  .+ Teuchos .+ Thyra .+ Stratimikos .+"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableThyraCoreLibs_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_ThyraCoreLibs:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Shards Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML RBGen 15"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_DefaultTPLs
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK DUMMY 3"
    "Final set of non-enabled TPLs:  MPI Boost Scotch METIS ParMETIS CppUnit ADOLC ADIC TVMET y12m SuperLUDist SuperLU UMFPACK AMD PETSC MUMPS 16"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoost
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoostParMETIS
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
    -DTPL_ENABLE_ParMETIS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost ParMETIS DUMMY 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS
  ARGS  -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting TPL_ENABLE_ParMETIS=ON because Zoltan_ENABLE_ParMETIS=ON"
    "Final set of enabled TPLs:  ParMETIS 1"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS_TplDisableParMETIS
  ARGS -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON -DTPL_ENABLE_ParMETIS=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Setting Zoltan_ENABLE_ParMETIS=OFF which was ON since TPL_ENABLE_ParMETIS=OFF"
    "Final set of enabled TPLs:  0"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnablePanzer_DefaultMPI
  ARGS  -DTrilinos_ENABLE_Panzer=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_Panzer=OFF which was 'ON' because Panzer has a required library dependence on disabled TPL MPI!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnablePanzer_EnableMPI
  ARGS  -DTrilinos_ENABLE_Panzer=ON
    -DTPL_ENABLE_MPI:BOOL=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  .*Panzer"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnablePanzer_DisableMPI
  ARGS  -DTrilinos_ENABLE_Panzer=ON
    -DTPL_ENABLE_MPI:BOOL=OFF
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "ERROR: Setting Trilinos_ENABLE_Panzer=OFF which was 'ON' because Panzer has a required library dependence on disabled TPL MPI!"
    "Configuring incomplete, errors occurred!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnablePanzer_DisableMPI_DisableEnabledDownstreamPackages
  ARGS  -DTrilinos_ENABLE_Panzer=ON
    -DTPL_ENABLE_MPI:BOOL=OFF
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
    -DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Setting Trilinos_ENABLE_Panzer=OFF which was 'ON' because Panzer has a required library dependence on disabled TPL MPI but Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON"
    "Final set of enabled packages:  0"
    "Final set of non-enabled packages:  .+ Panzer .+"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTpetra_EnableAllForwardPackages
  ARGS  -DTrilinos_ENABLE_Tpetra=ON
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra 8"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt ThyraTpetra Thyra 12"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryTestedCode_EnableAllOptionalPackges
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON -DTrilinos_ENABLE_SECONDARY_TESTED_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx 22"
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


# Deprecated ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE!
CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryStableCode_EnableAllOptionalPackges_BackwardCompatible
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Set var Trilinos_ENABLE_SECONDARY_STABLE_CODE is deprecated!  Use Trilinos_ENABLE_SECONDARY_TESTED_CODE instead!"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx 22"
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


# Test that deprecated Trilinos_ENABLE_SECONDARY_STABLE_CODE is not set when
# Trilinos_ENABLE_SECONDARY_TESTED_CODE is also being set.
CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryTestedCode_EnableSecondaryStableCode_EnableAllOptionalPackges
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON  -DTrilinos_ENABLE_SECONDARY_TESTED_CODE=ON
    -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Set deprecated Trilinos_ENABLE_SECONDARY_STABLE_CODE and"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem1
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem1
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack ML 3"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Setting Stratimikos_ENABLE_ML=OFF because Stratimikos has an optional library dependence on disabled package ML"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra AztecOO Galeri Amesos Stratimikos 11"
    "Final set of non-enabled packages:  TrilinosFramework Shards Tpetra Stokhos Sacado Isorropia Intrepid Ifpack ML Belos RBGen Phalanx Panzer 13"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack 2"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 13"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 16"
    "Final set of non-enabled packages:  TrilinosFramework Shards Tpetra Stokhos Sacado Intrepid Ifpack Belos RBGen Phalanx Panzer 11"
    "Final set of non-enabled SE packages:  TrilinosFramework Shards Tpetra Stokhos Sacado ThyraGoodStuff ThyraCrazyStuff ThyraTpetra Intrepid Ifpack Belos RBGen Phalanx Panzer 14"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableIfpack_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Ifpack=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: User has set Trilinos_ENABLE_Ifpack=ON but the" 
    "package Ifpack is not supported on this platform type 'BadSystem2'"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 14"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 17"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableTeko
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because Teko has a required dependence on Teuchos"
    "Setting Trilinos_ENABLE_Epetra=ON because Teko has a required dependence on Epetra"
    "Setting Trilinos_ENABLE_Thyra=ON because Teko has a required dependence on Thyra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 16"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Panzer Stalix 10"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL 3"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableTeko_InconsistentExtraReposFile
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_EXTRAREPOS_FILE=${CMAKE_CURRENT_SOURCE_DIR}/NoPackagesExtraReposFile.cmake
    -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=ON
    #-DTrilinos_VERBOSE_CONFIGURE=ON
    #-DTRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG=ON
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
    "Explicitly enabled packages on input .* Teko 1"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DefaultExtraReposFile_EnableAllPackages_TraceFileProcessing
  ARGS 
    -DTrilinos_EXTRAREPOS_FILE=${MOCK_TRILINOS_DIR}/cmake/ExtraRepositoriesList.cmake
    -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=TRUE
    -DTrilinos_TRACE_FILE_PROCESSING=ON
    -DTrilinos_ENABLE_ALL_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "-- File Trace: PROJECT    INCLUDE    .*/MockTrilinos/Version.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/cmake/CallbackSetupExtraOptions.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/PackagesList.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/TPLsList.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/preCopyrightTrilinos/PackagesList.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/preCopyrightTrilinos/TPLsList.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/extraTrilinosRepo/PackagesList.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/extraTrilinosRepo/TPLsList.cmake"
    "-- File Trace: REPOSITORY INCLUDE    .*/MockTrilinos/cmake/RepositoryDependenciesSetup.cmake"
    "-- File Trace: PROJECT    INCLUDE    .*/MockTrilinos/cmake/ProjectDependenciesSetup.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/cmake/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/teuchos/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/src/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/src/good_stuff/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/src/crazy_stuff/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/adapters/epetra/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/adapters/epetraext/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/packages/thyra/adapters/tpetra/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/preCopyrightTrilinos/teko/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/preCopyrightTrilinos/stalix/cmake/Dependencies.cmake"
    "-- File Trace: PACKAGE    INCLUDE    .*/MockTrilinos/extraTrilinosRepo/extrapack/cmake/Dependencies.cmake"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableStalix
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_ENABLE_Stalix=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Stalix 1"
    "Explicitly disabled packages on input .* Stokhos 1"
    "WARNING: Not Setting Trilinos_ENABLE_ThyraGoodStuff=ON even though Thyra has an optional dependence on ThyraGoodStuff because Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF"
    "Setting Trilinos_ENABLE_Tpetra=ON because Stalix has a required dependence on Tpetra"
    "Setting Trilinos_ENABLE_Thyra=ON because Stalix has a required dependence on Thyra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Sacado Thyra Stalix 10"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 16"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidPackageNameError
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DSHOW_INVALID_PACKAGE_NAME_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'Anasazi' is listed as a dependency of the package"
    "'Teko' is in the list 'Trilinos_SE_PACKAGES' but the package 'Anasazi' is"
    "either not defined or is listed later in the package order.  Check the"
    "pelling of 'Anasazi' or see how it is listed in"
    "Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS in relationship to 'Teko'"
    "Configuring incomplete, errors occurred"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidPackageSelfDependency
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DSHOW_SELF_PACKAGE_DEPENDENCY_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'Teko' is listed as a dependency itself in the list"
    " 'LIB_REQUIRED_DEP_PACKAGES'!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidSubpackageSelfDependency
  ARGS
    -DSHOW_SELF_SUBPACKAGE_DEPENDENCY_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'ThyraTpetra' is listed as a dependency itself in the"
     " list 'LIB_REQUIRED_DEP_PACKAGES'!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  MissingExtraRepositoryFail
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=MissingExtraRepo
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "ERROR: Skipping missing extra repo 'MissingExtraRepo' packages list file"
    "ERROR: Skipping missing extra repo 'MissingExtraRepo' TPLs list file"
    "Configuring incomplete, errors occurred"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  MissingExtraRepositoryIgnore
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=MissingExtraRepo
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=ON
    -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=
    -DTRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG=TRUE
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING!  Ignoring missing extra repo 'MissingExtraRepo' packages list file ..+/package_arch/UnitTests/MockTrilinos/MissingExtraRepo/PackagesList.cmake. on request!"
    "WARNING!  Ignoring missing extra repo 'MissingExtraRepo' TPLs list file ..+/package_arch/UnitTests/MockTrilinos/MissingExtraRepo/TPLsList.cmake. on request!"
    "Configuring done"
    "Generating done"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  NativeRepositoriesList1_EnableAllPackages_ST
  ARGS -DTrilinos_NATIVE_REPO_FILE_OVERRRIDE=cmake/NativeRepositoriesList1.cmake
    -DTrilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=ON
    -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Reading the list of packages from .*/MockTrilinos/PackagesList.cmake"
    "Reading the list of packages from .*/MockTrilinos/preCopyrightTrilinos/PackagesList.cmake"
    "Reading the list of packages from .*/MockTrilinos/extraTrilinosRepo/PackagesList.cmake"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Teko ExtraPack 24"
    "Configuring done"
    "Generating done"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_EnableExtraPack
  ARGS "-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos,extraTrilinosRepo"
    -DTrilinos_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because ExtraPack has a required dependence on Teuchos"
    "Setting Trilinos_ENABLE_Teko=ON because ExtraPack has a required dependence on Teko"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even though ExtraPack has an optional dependence on Sacado because Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko ExtraPack 17"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Panzer Stalix 10"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL2 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_EnableExtraPack_EnableExtraTPL1
  ARGS "-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos,extraTrilinosRepo"
    -DTrilinos_ENABLE_ExtraPack=ON -DTPL_ENABLE_ExtraTPL1=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL1 ExtraTPL2 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_ReplacePackagesList
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DpreCopyrightTrilinos_PACKAGES_LIST_FILE=cmake/Reduced_preCopyrightTrilinos_PackagesList.cmake
  PASS_REGULAR_EXPRESSION_ALL
    "Reading a list of extra packages from .*/MockTrilinos/cmake/Reduced_preCopyrightTrilinos_PackagesList.cmake"
    "Trilinos_NUM_PACKAGES='24'"
    "Final set of non-enabled SE packages:  .* Phalanx Panzer Stalix 31"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_CI_EnableTeko
  ARGS -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Continuous
   -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=ON
   -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because Teko has a required dependence on Teuchos"
    "Setting Trilinos_ENABLE_Epetra=ON because Teko has a required dependence on Epetra"
    "Setting Trilinos_ENABLE_Thyra=ON because Teko has a required dependence on Thyra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 16"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Panzer Stalix 10"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_Nightly_EnableExtraPack
  ARGS  -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Nightly
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=ON
    -DTrilinos_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because ExtraPack has a required dependence on Teuchos"
    "Setting Trilinos_ENABLE_Teko=ON because ExtraPack has a required dependence on Teko"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko ExtraPack 17"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Panzer Stalix 10"
  )
