# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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
INCLUDE(AdvancedSet)


ADVANCED_SET(Trilinos_CTEST_COMMAND ctest CACHE STRING
  "CTest executable used to drive unit tests")

# The mock test project 
SET(MOCK_PROJECT_NAME Trilinos)

SET(MOCK_TRILINOS_DIR ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_MOCK_TRILINOS_DIR})


FUNCTION(CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE TEST_NAME)

  #MESSAGE("CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE: ${TEST_CASE} ${ARGN}")

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "COPY_FILES;ARGS;PASS_REGULAR_EXPRESSION_ALL;COMPARE_DEPS_FILE"
     #options
     "SKIP_STANDARD_DEPS_CHECK"
     ${ARGN}
     )
  
  SET(TEST_CORE_NAME CTestDriver_${TEST_NAME})
  SET(TEST_FULL_NAME ${PACKAGE_NAME}_${TEST_CORE_NAME})

  SET(TEST_BLOCKS "")
  SET(TEST_IDX 0)

  IF (PARSE_COPY_FILES)

    LIST(LENGTH  PARSE_COPY_FILES  COPY_FILES_LEN)
    MATH(EXPR NUM_FILE_COPIES "${COPY_FILES_LEN}/2")
    MATH(EXPR COPY_FILES_LAST_IDX "${NUM_FILE_COPIES}-1")

    SET(NUM_COPY_FILE_ARGS 2)
    SET(FILE_FROM_OFFSET 0)
    SET(FILE_TO_OFFSET 1)

    FOREACH(COPY_FILE_IDX RANGE ${COPY_FILES_LAST_IDX})

      # FILE_FROM
      MATH(EXPR FILE_FROM_IDX "${COPY_FILE_IDX}*${NUM_COPY_FILE_ARGS}+${FILE_FROM_OFFSET}")
      LIST(GET PARSE_COPY_FILES ${FILE_FROM_IDX} FILE_FROM )
      #PRINT_VAR(FILE_FROM)
 
      # FILE_TO
      MATH(EXPR FILE_TO_IDX "${COPY_FILE_IDX}*${NUM_COPY_FILE_ARGS}+${FILE_TO_OFFSET}")
      LIST(GET PARSE_COPY_FILES ${FILE_TO_IDX} FILE_TO )
      #PRINT_VAR(FILE_TO)

      SET(TEST_BLOCKS ${TEST_BLOCKS}
        TEST_${TEST_IDX} CMND ${CMAKE_COMMAND}
          ARGS -E copy ${FILE_FROM} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_FULL_NAME}/${FILE_TO}
        )

      MATH(EXPR TEST_IDX "${TEST_IDX}+1")

    ENDFOREACH()

  ENDIF()

  SET(TEST_BLOCKS ${TEST_BLOCKS}
    TEST_${TEST_IDX} CMND env
      ARGS
        ${MOCK_PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        CTEST_DEPENDENCY_HANDLING_UNIT_TESTING=ON
        ${MOCK_PROJECT_NAME}_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}
        ${MOCK_PROJECT_NAME}_PACKAGES=  # Set back to default
        CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES= # Set back to default
        CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=OFF # Set back to default
        CTEST_DO_UPDATES=ON # Set back to default
        CTEST_DASHBOARD_ROOT=  # Must be empty
        CTEST_DO_SUBMIT=OFF # Must be off
        CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=OFF  # Must be off!
        ${PARSE_ARGS}
        ${CMAKE_CTEST_COMMAND} -S ${CMAKE_CURRENT_SOURCE_DIR}/unit_test_driver.cmake
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    )

  IF (NOT PARSE_SKIP_STANDARD_DEPS_CHECK)
    SET(TEST_BLOCKS ${TEST_BLOCKS}
      "Dumping the XML dependencies file .*/${TEST_FULL_NAME}/TrilinosPackageDependencies.xml"
      "CDASH_SUBPROJECT_XML_FILE='.*/${TEST_FULL_NAME}/CDashSubprojectDependencies.xml'"
      )
    IF (PYTHON_EXECUTABLE)
      SET(TEST_BLOCKS ${TEST_BLOCKS}
        "Dumping the CDash XML dependencies file .*/${TEST_FULL_NAME}/CDashSubprojectDependencies.xml"
        )
    ENDIF()
  ENDIF()

  MATH(EXPR TEST_IDX "${TEST_IDX}+1")

  IF (PARSE_COMPARE_DEPS_FILE)
    SET(TEST_BLOCKS ${TEST_BLOCKS}
      TEST_${TEST_IDX} CMND ${CMAKE_COMMAND}
        ARGS -E compare_files
          "${PARSE_COMPARE_DEPS_FILE}"
          "${CMAKE_CURRENT_BINARY_DIR}/${TEST_FULL_NAME}/TrilinosPackageDependencies.xml"
      )
  ENDIF()

  TRIBITS_ADD_ADVANCED_TEST( ${TEST_CORE_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    ${TEST_BLOCKS}
    )

ENDFUNCTION()


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllDefault
  PASS_REGULAR_EXPRESSION_ALL
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "Final set of non-enabled packages:  Stokhos Sacado Belos Phalanx Panzer 5"
  COMPARE_DEPS_FILE
     "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependencies.Trilinos.xml"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllSecondaryStable
  ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 23"
    "Final set of non-enabled packages:  Stokhos 1"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTeuchos
  ARGS Trilinos_PACKAGES=Teuchos Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Experimental
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos 1"
    "Final set of non-enabled packages:  TrilinosFramework RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Stokhos Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 23"
    "Processing current package Teuchos: libs='ON', tests='ON'"
    "CONFIGURE_OPTIONS = '-DTrilinos_TRIBITS_DIR=.+-DCTEST_USE_LAUNCHERS:BOOL=1.-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON.-DTrilinos_ENABLE_TESTS:BOOL=ON.-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=-DummyErrFlags.-DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON.-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON.-DMPI_EXEC_MAX_NUMPROCS:STRING=4.-DTrilinos_EXTRAREPOS_FILE:STRING=.*/MockTrilinos/cmake.ExtraRepositoriesList.cmake.-DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES:BOOL=ON.-DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Experimental.-DTrilinos_ENABLE_Teuchos:BOOL=ON'"
  )
# NOTE: The above test pins down the form of the configure options passed to
# different package configures to make sure everything in there.


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra
  ARGS Trilinos_PACKAGES=Thyra Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Experimental
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra 8"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 16"
    "Processing current package Teuchos: libs='ON', tests=''"
    "CONFIGURE_OPTIONS = '-DTrilinos_TRIBITS_DIR=.+-DCTEST_USE_LAUNCHERS:BOOL=1.-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON.-DTrilinos_ENABLE_TESTS:BOOL=.-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=-DummyErrFlags.-DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON.-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON.-DMPI_EXEC_MAX_NUMPROCS:STRING=4.-DTrilinos_EXTRAREPOS_FILE:STRING=.+/tribits/package_arch/UnitTests/MockTrilinos/cmake/ExtraRepositoriesList.cmake.-DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES:BOOL=ON.-DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Experimental.-DTrilinos_ENABLE_Teuchos:BOOL=ON'"
    "Processing current package Epetra: libs='ON', tests=''"
    "Processing current package Thyra: libs='ON', tests='ON'"
    "CONFIGURE_OPTIONS = '-DTrilinos_TRIBITS_DIR=.+-DCTEST_USE_LAUNCHERS:BOOL=1.-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON.-DTrilinos_ENABLE_TESTS:BOOL=ON.-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=-DummyErrFlags.-DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON.-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON.-DMPI_EXEC_MAX_NUMPROCS:STRING=4.-DTrilinos_EXTRAREPOS_FILE:STRING=.+/tribits/package_arch/UnitTests/MockTrilinos/cmake/ExtraRepositoriesList.cmake.-DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES:BOOL=ON.-DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Experimental.-DTrilinos_ENABLE_EpetraExt:BOOL=.-DTrilinos_ENABLE_Thyra:BOOL=ON'"
    "Processing current package Stokhos: libs='OFF', tests=''"
    "Processing current package Isorropia: libs='', tests=''"
  )
# NOTE: The above test pins down the form of the configure options passed to
# different package configures to make sure everything in there.


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_ExcludeEpetra
  ARGS Trilinos_PACKAGES=Thyra Trilinos_EXCLUDE_PACKAGES=Epetra
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Thyra 1"
    "Explicitly disabled packages on input .*  Epetra Stokhos 2"
    "Final set of enabled packages:  Teuchos RTOp Tpetra Thyra 4"
    "Final set of non-enabled packages:  TrilinosFramework Epetra Zoltan Shards Triutils EpetraExt Stokhos Sacado Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 20"
    "Processing current package Teuchos: libs='ON', tests=''"
    "Processing current package Epetra: libs='OFF', tests=''"
    "Processing current package Zoltan: libs='', tests=''"
    "Processing current package EpetraExt: libs='OFF', tests=''"
    "Processing current package Thyra: libs='ON', tests='ON'"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAll_EnablePhalanx
  ARGS Trilinos_ADDITIONAL_PACKAGES=Phalanx
  PASS_REGULAR_EXPRESSION_ALL
    "Enabling explicitly set package Phalanx"
    "Explicitly enabled packages on input .*  Phalanx 1"
    "Explicitly disabled packages on input .*  Stokhos 1"
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen Phalanx 21"
    "Final set of non-enabled packages:  Stokhos Belos Panzer 3"
    "Processing current package Stokhos: libs='OFF', tests=''"
    "Processing current package RBGen: libs='ON', tests='ON'"
    "Processing current package Phalanx: libs='ON', tests='ON'"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAll_EnableStokhos
  ARGS Trilinos_ADDITIONAL_PACKAGES=Stokhos
  PASS_REGULAR_EXPRESSION_ALL
    "Enabling explicitly set package Stokhos"
    "Explicitly enabled packages on input .*  Stokhos 1"
    "Explicitly disabled packages on input .*  0"
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Stokhos Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 20"
    "Final set of non-enabled packages:  Sacado Belos Phalanx Panzer 4"
    "Processing current package Stokhos: libs='ON', tests='ON'"
    "Processing current package RBGen: libs='ON', tests='ON'"
    "Processing current package Phalanx: libs='', tests=''"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllSecondaryStable_EnableStokhos
  ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON Trilinos_ADDITIONAL_PACKAGES=Stokhos
  PASS_REGULAR_EXPRESSION_ALL
    "Enabling explicitly set package Stokhos"
    "Explicitly enabled packages on input .*  Stokhos 1"
    "Explicitly disabled packages on input .*  0"
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Stokhos Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer 24"
    "Final set of non-enabled packages:  0"
    "Processing current package Stokhos: libs='ON', tests='ON'"
    "Processing current package RBGen: libs='ON', tests='ON'"
    "Processing current package Phalanx: libs='ON', tests='ON'"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  CI_EnableAllSecondaryStable
  ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Continuous
  PASS_REGULAR_EXPRESSION_ALL
    "EXECUTE_PROCESS.+git.+preCopyrightTrilinos.+preCopyrightTrilinos.update.out"
    "Explicitly enabled packages on input .*  0"
    "Explicitly disabled packages on input .*  Stokhos Stalix 2"
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 24"
    "Final set of non-enabled packages:  Stokhos Stalix 2"
    "CONFIGURE_OPTIONS.+-DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Continuous.+"
    "Processing current package Teuchos: libs='ON', tests='ON'"
    "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Teuchos:BOOL=ON"
    "Processing current package Phalanx: libs='ON', tests='ON'"
    "Processing current package Stalix: libs='OFF', tests=''"
    "Package Stalix is disabled, skipping configure, build, test"
  COMPARE_DEPS_FILE
     "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependencies.Continuos.xml"
  )


# Check with an extra repo enabled and assert all of the internal logic
CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  CI_EnableAllSecondaryStable_EnableStalix
  ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Continuous
    Trilinos_ADDITIONAL_PACKAGES=Stalix
  PASS_REGULAR_EXPRESSION_ALL
    "preCopyrightTrilinos: Doing GIT update from URL 'url1:/git/preCopyrightTrilinos' to dir '.*/MockTrilinos/preCopyrightTrilinos'"
    "EXECUTE_PROCESS.COMMAND./somebasedir/git.reset.--hard.WORKING_DIRECTORY..*/MockTrilinos/preCopyrightTrilinos.OUTPUT_FILE..*/tribits/ctest/UnitTests/TrilinosFramework_CTestDriver_CI_EnableAllSecondaryStable_EnableStalix/preCopyrightTrilinos.reset.out."
    "EXECUTE_PROCESS.COMMAND./somebasedir/git.clean.-df.WORKING_DIRECTORY..*/MockTrilinos/preCopyrightTrilinos.OUTPUT_FILE..*/tribits/ctest/UnitTests/TrilinosFramework_CTestDriver_CI_EnableAllSecondaryStable_EnableStalix/preCopyrightTrilinos.clean.out."
    "EXECUTE_PROCESS.COMMAND./somebasedir/git.pull.WORKING_DIRECTORY..*/MockTrilinos/preCopyrightTrilinos.OUTPUT_FILE..*/cmake/tribits/ctest/UnitTests/TrilinosFramework_CTestDriver_CI_EnableAllSecondaryStable_EnableStalix/preCopyrightTrilinos.update.out."
    "EXECUTE_PROCESS.COMMAND;/somebasedir/git.log.--pretty=format:.*.-C.HEAD..ORIG_HEAD.OUTPUT_VARIABLE.CLONE_OR_PULL_OUTPUT.WORKING_DIRECTORY;.*/tribits/package_arch/UnitTests/MockTrilinos/preCopyrightTrilinos."
    "Dakota: Doing initial SVN clone/checkout from URL 'url3/somedirpath/trunk' to dir 'packages/TriKota/Dakota'"
    "EXECUTE_PROCESS.COMMAND./someotherbasedir/svn.checkout.url3/somedirpath/trunk.packages/TriKota/Dakota.WORKING_DIRECTORY..*/MockTrilinos.OUTPUT_FILE..*/cmake/tribits/ctest/UnitTests/TrilinosFramework_CTestDriver_CI_EnableAllSecondaryStable_EnableStalix/Dakota.update.out."
    "EXECUTE_PROCESS.echo.OUTPUT_VARIABLE.CLONE_OR_PULL_OUTPUT.WORKING_DIRECTORY.*/tribits/package_arch/UnitTests/MockTrilinos/packages/TriKota/Dakota."
    ".... preCopyrightTrilinos updates."
    "... Dakota cloned version."
    "Adding extra Continuous repository preCopyrightTrilinos"
    "Adding extra Continuous repository Dakota"
    "Trilinos_NUM_PACKAGES='24'"
    "Trilinos_NUM_TPLS='19'"
    "Reading a list of extra packages from .*/MockTrilinos/preCopyrightTrilinos/PackagesList.cmake"
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
    "Skipping reading packages and TPLs for extra repo Dakota because marked NOPACKAGES"
    "Explicitly enabled packages on input .*  Stalix 1"
    "Explicitly disabled packages on input .*  Stokhos 1"
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko Stalix 25"
    "Final set of non-enabled packages:  Stokhos 1"
    "CTEST_NOTES_FILES_WO_CACHE=.*/Updates.txt"
    "Processing current package Phalanx: libs='ON', tests='ON'"
    "Processing current package Stalix: libs='ON', tests='ON'"
  )


CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
  Nightly_EnableAllSecondaryStable
  ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly
  PASS_REGULAR_EXPRESSION_ALL
    "EXECUTE_PROCESS.+git.+preCopyrightTrilinos.+preCopyrightTrilinos.update.out"
    "EXECUTE_PROCESS.+git.+extraTrilinosRepo.+extraTrilinosRepo.update.out"
    "Explicitly enabled packages on input .*  0"
    "Explicitly disabled packages on input .*  Stokhos Stalix 2"
    "Enabling all packages that are not currently disabled because of Trilinos_ENABLE_ALL_PACKAGES=ON"
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko ExtraPack 25"
    "Final set of non-enabled packages:  Stokhos Stalix 2"
    "Trilinos_NUM_ENABLED_PACKAGES=25:  Configuring packages!"
    "Processing current package ExtraPack: libs='ON', tests='ON'"
    "Configuring TRIBITS_PACKAGE='ExtraPack'"
  )


#
# Test enable logic based on modified files/packages
#


IF (PYTHON_EXECUTABLE)

  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedNone
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
    PASS_REGULAR_EXPRESSION_ALL
      "Final set of enabled packages:  0"
      "MODIFIED_PACKAGES_LIST='':  No modified packages to justify continuous integration test iteration"
    SKIP_STANDARD_DEPS_CHECK
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedTeuchos
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Teuchos.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
    ARGS  Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous  CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_UPDATE_RETURN_VAL=1
    PASS_REGULAR_EXPRESSION_ALL
      "Directly modified or failing non-disabled packages that need to be tested:  Teuchos 1"
      "Explicitly enabled packages on input .*  Teuchos 1"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 23"
      "Final set of non-enabled packages:  TrilinosFramework Stokhos Stalix 3"
      "CTEST_UPDATE.+ returned '1'"
      "Processing current package Teuchos: libs='ON', tests='ON'"
      "Processing current package Teko: libs='ON', tests='ON'"
      "Processing current package Stalix: libs='OFF', tests=''"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_Nightly_EnableOnlyModified_ModifiedTeuchos
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Teuchos.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.extraTrilinosRepo.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE
      Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly
      CTEST_UPDATE_RETURN_VAL=1
    PASS_REGULAR_EXPRESSION_ALL
      "Directly modified or failing non-disabled packages that need to be tested:  Teuchos 1"
      "Explicitly enabled packages on input .*  Teuchos 1"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko ExtraPack 24"
      "Final set of non-enabled packages:  TrilinosFramework Stokhos Stalix 3"
      "CTEST_UPDATE.+ returned '1'"
      "Processing current package Teuchos: libs='ON', tests='ON'"
      "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Teuchos:BOOL=ON"
      "Processing current package Zoltan: libs='ON', tests=''"
      "Processing current package Teko: libs='ON', tests='ON'"
      "Processing current package Stalix: libs='OFF', tests=''"
      "Processing current package ExtraPack: libs='ON', tests='ON'"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedThyra_NoEnableImplictlyEnabled
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Thyra.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.extraTrilinosRepo.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
      CTEST_UPDATE_RETURN_VAL=1
    PASS_REGULAR_EXPRESSION_ALL
      "Directly modified or failing non-disabled packages that need to be tested:  Thyra 1"
      "Explicitly enabled packages on input .*  Thyra 1"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos Phalanx Panzer Teko 22"
      "Final set of non-enabled packages:  TrilinosFramework Stokhos RBGen Stalix 4"
      "CTEST_UPDATE.+ returned '1'"
      "Processing current package Teuchos: libs='ON', tests=''"
      "Not enabling implicitly enabled package Teuchos on request"
      "Processing current package Zoltan: libs='ON', tests=''"
      "Not enabling implicitly enabled package Zoltan on request"
      "Processing current package Thyra: libs='ON', tests='ON'"
      "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Thyra:BOOL=ON"
      "Processing current package Teko: libs='ON', tests='ON'"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedThyraCoreLibs_NoEnableImplictlyEnabled
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.ThyraCoreLibs.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.extraTrilinosRepo.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
      CTEST_UPDATE_RETURN_VAL=1
    PASS_REGULAR_EXPRESSION_ALL
      "Directly modified or failing non-disabled packages that need to be tested:  ThyraCoreLibs 1"
      "Explicitly enabled packages on input .*  0"
      "Explicitly enabled SE packages on input .*  ThyraCoreLibs 1"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos Phalanx Panzer Teko 22"
      "Final set of non-enabled packages:  TrilinosFramework Stokhos RBGen Stalix 4"
      "CTEST_UPDATE.+ returned '1'"
      "Processing current package Teuchos: libs='ON', tests=''"
      "Not enabling implicitly enabled package Teuchos on request"
      "Processing current package Zoltan: libs='ON', tests=''"
      "Not enabling implicitly enabled package Zoltan on request"
      "Processing current package Thyra: libs='ON', tests='ON'"
      "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Thyra:BOOL=ON"
      "Processing current package Teko: libs='ON', tests='ON'"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedTeko
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Teko.txt modifiedFiles.preCopyrightTrilinos.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE
      CTEST_UPDATE_RETURN_VAL=0 # Main update returns nothing for the main Trilinos repo!
    PASS_REGULAR_EXPRESSION_ALL
      "Directly modified or failing non-disabled packages that need to be tested:  Teko 1"
      "Explicitly enabled packages on input .*  Teko 1"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Belos Stratimikos Teko 17"
      "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid RBGen Phalanx Panzer Stalix 9"
      "CTEST_UPDATE.+ returned '0'"
      "Processing current package Teuchos: libs='ON', tests=''"
      "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Teuchos:BOOL=ON"
      "Processing current package Teko: libs='ON', tests='ON'"
      "Processing current package Stalix: libs='OFF', tests=''"
    )
  # NOTE: The above unit test is critical because it shows that Teko will get
  # enabled and processed even though the main git update of the main Trilinos
  # repo returns no updated files.  This ensures that if only extra repo
  # packages get modified then they will get picked up in CI testing and
  # processed correctly.


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedTeuchosTeko
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Teuchos.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Teko.txt modifiedFiles.preCopyrightTrilinos.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_UPDATE_RETURN_VAL=2
    PASS_REGULAR_EXPRESSION_ALL
      "MODIFIED_PACKAGES_LIST='Teuchos.Teko'"
      "Enabling modified package: Teuchos"
      "Enabling modified package: Teko"
      "Directly modified or failing non-disabled packages that need to be tested:  Teuchos Teko 2"
      "Explicitly enabled packages on input .*  Teuchos Teko 2"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 23"
      "Final set of non-enabled packages:  TrilinosFramework Stokhos Stalix 3"
      "CTEST_UPDATE.+ returned '2'"
      "Processing current package Teuchos: libs='ON', tests='ON'"
      "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Teuchos:BOOL=ON"
      "Processing current package Teko: libs='ON', tests='ON'"
      "CONFIGURE_OPTIONS.+Trilinos_ENABLE_Teko:BOOL=ON"
      "Processing current package Stalix: libs='OFF', tests=''"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedStokhosStalix_FailedNone
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Stokhos.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Stalix.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/failedPackages.none.txt failedPackages.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_UPDATE_RETURN_VAL=1
    PASS_REGULAR_EXPRESSION_ALL
      "MODIFIED_PACKAGES_LIST='Stokhos.Stalix'"
      "FAILING_PACKAGES_LIST=''"
      "Not enabling explicitly disabled modified package: Stokhos"
      "Not enabling explicitly disabled modified package: Stalix"
      "Directly modified or failing non-disabled packages that need to be tested:  0"
      "Explicitly enabled packages on input .*  0"
      "Explicitly disabled packages on input .*  Stokhos Stalix 2"
      "Final set of enabled packages:  0"
      "Final set of non-enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Stokhos Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko Stalix 26"
      "Trilinos_NUM_ENABLED_PACKAGES=0:  Exiting the script!"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedTeko_FailedThyraRBGen
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.Teko.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/failedPackages.Thyra.RBGen.txt failedPackages.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_UPDATE_RETURN_VAL=0 # Main update returns nothing for the main Trilinos repo!
    PASS_REGULAR_EXPRESSION_ALL
      "MODIFIED_PACKAGES_LIST='Teko'"
      "FAILING_PACKAGES_LIST='Thyra.RBGen'"
      "Enabling modified package: Teko"
      "Enabling previously failing package: Thyra"
      "Enabling previously failing package: RBGen"
      "Directly modified or failing non-disabled packages that need to be tested:  Thyra RBGen Teko 3"
      "Explicitly enabled packages on input .*  Thyra RBGen Teko 3"
      "Explicitly disabled packages on input .*  Stokhos Stalix 2"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 23"
      "Final set of non-enabled packages:  TrilinosFramework Stokhos Stalix 3"
      "CTEST_UPDATE.+ returned '0'"
      "Processing current package Teuchos: libs='ON', tests=''"
      "Processing current package Epetra: libs='ON', tests=''"
      "Processing current package Thyra: libs='ON', tests='ON'"
      "Processing current package RBGen: libs='ON', tests='ON'"
      "Processing current package Stokhos: libs='OFF', tests=''"
      "Processing current package Teko: libs='ON', tests='ON'"
      "Processing current package Stalix: libs='OFF', tests=''"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedNone_FailedThyraRBGen
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/failedPackages.Thyra.RBGen.txt failedPackages.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_UPDATE_RETURN_VAL=0 # Main update returns nothing for the main Trilinos repo!
    PASS_REGULAR_EXPRESSION_ALL
      "CTEST_UPDATE.+ returned '0'"
      "MODIFIED_PACKAGES_LIST=''"
      "FAILING_PACKAGES_LIST='Thyra.RBGen'"
      "Directly modified or failing non-disabled packages that need to be tested:  Thyra RBGen 2"
      "Explicitly enabled packages on input .*  Thyra RBGen 2"
      "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 23"
      "MODIFIED_PACKAGES_LIST='':  No modified packages to justify continuous integration test iteration"
    )


  CREATE_CTEST_DEPENDENCY_HANDLING_TEST_CASE(
    CI_EnableOnlyModified_ModifiedNone_FailedStokhosStalix
    COPY_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/modifiedFiles.empty.txt modifiedFiles.preCopyrightTrilinos.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/failedPackages.Stokhos.Stalix.txt failedPackages.txt
    ARGS Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON
      CTEST_TEST_TYPE=Continuous CTEST_ENABLE_MODIFIED_PACKAGES_ONLY=ON
      CTEST_UPDATE_RETURN_VAL=0 # Main update returns nothing for the main Trilinos repo!
    PASS_REGULAR_EXPRESSION_ALL
      "CTEST_UPDATE.+ returned '0'"
      "MODIFIED_PACKAGES_LIST=''"
      "FAILING_PACKAGES_LIST='Stokhos.Stalix'"
      "Not enabling explicitly disabled previously failing package: Stokhos"
      "Not enabling explicitly disabled previously failing package: Stalix"
      "Directly modified or failing non-disabled packages that need to be tested:  0"
      "Explicitly enabled packages on input .*  0"
      "Explicitly disabled packages on input .*  Stokhos Stalix 2"
      "Final set of enabled packages:  0"
      "Final set of non-enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Stokhos Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko Stalix 26"
      "MODIFIED_PACKAGES_LIST='':  No modified packages to justify continuous integration test iteration"
    )

ENDIF(PYTHON_EXECUTABLE)

# ToDo: Test ???
