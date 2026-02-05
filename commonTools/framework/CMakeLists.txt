TRIBITS_PACKAGE(TrilinosFrameworkTests)


add_subdirectory(github_issue_creator)


TRIBITS_ADD_ADVANCED_TEST( ProjectCiFileChangeLogic_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND ${Python3_EXECUTABLE}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ProjectCiFileChangeLogic_UnitTests.py -v
    PASS_REGULAR_EXPRESSION "OK"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  )


TRIBITS_ADD_ADVANCED_TEST( get-changed-trilinos-packages-helpers-unit-tests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND
      "${CMAKE_CURRENT_SOURCE_DIR}/get-changed-trilinos-packages-helpers-unit-tests.sh"
    PASS_REGULAR_EXPRESSION "OK"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  )


TRIBITS_ADD_ADVANCED_TEST( clean_all_jobs_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  EXCLUDE_IF_NOT_TRUE TFW_Python_Testing
  TEST_0 CMND ${Python3_EXECUTABLE}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/clean_workspace/unittests/test_clean_all_jobs.py -v
    PASS_REGULAR_EXPRESSION "OK"
  )


TRIBITS_ADD_ADVANCED_TEST( clean_sentinel_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  EXCLUDE_IF_NOT_TRUE TFW_Python_Testing
  TEST_0 CMND ${Python3_EXECUTABLE}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/clean_workspace/unittests/test_clean_sentinel.py -v
    PASS_REGULAR_EXPRESSION "OK"
  )


TRIBITS_ADD_ADVANCED_TEST( clean_workspace_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  EXCLUDE_IF_NOT_TRUE TFW_Python_Testing
  TEST_0 CMND ${Python3_EXECUTABLE}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/clean_workspace/unittests/test_clean_workspace.py -v
    PASS_REGULAR_EXPRESSION "OK"
  )

# Run the discoverable unit tests under the packages/framework/pr_tools directory tree
TRIBITS_ADD_ADVANCED_TEST( PullRequestLinuxDriver_UnitTests
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    EXCLUDE_IF_NOT_TRUE TFW_Python_Testing
    TEST_0 CMND bash
        ARGS pr_python_test_driver.sh ${Python3_EXECUTABLE} ${PYTHON_PIP_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../packages/framework/pr_tools
        SKIP_CLEAN_WORKING_DIRECTORY
        PASS_REGULAR_EXPRESSION "=+ [0-9]+ passed"
        #ALWAYS_FAIL_ON_NONZERO_RETURN
)
# Note: pytest output looks like this:
# Passing tests:
# ========...======= 36 passed in 0.79s ====...======
# Failing tests:
# ========...======= 1 failed, 35 passed in 0.80s ======...=======



FUNCTION(create_get_changed_trilinos_packages_test  TEST_POSTFIX
  FILES_CHANGED
  CHANGED_PACKAGES_FULL_LIST
  CHANGED_PACKAGES_ST_LIST
  EXPECTED_ENABLES_REGEX
  EXPECTED_CTEST_LABELS_FOR_SUBPROJECTS_REGEX
  )

  SET(TEST_NAME "get_changed_trilinos_packages_test_${TEST_POSTFIX}")
  SET(FULL_TEST_NAME "${PACKAGE_NAME}_${TEST_NAME}")
  SET(ABS_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/${FULL_TEST_NAME}")

  # We can't use real git in these automated tests so we use the handy
  # mockprogram.py script instead of git.
  SET(GIT_MOCKPROGRAM_TXT_FILE
    "${CMAKE_CURRENT_BINARY_DIR}/${FULL_TEST_NAME}_git_mockprogram.txt")

  # We need to write the expected input and the output for
  FILE(WRITE "${GIT_MOCKPROGRAM_TXT_FILE}"
"MOCK_PROGRAM_INPUT: diff --name-only sha1-from..sha1-to
MOCK_PROGRAM_RETURN: 0
MOCK_PROGRAM_OUTPUT: ${FILES_CHANGED}
"
  )

  TRIBITS_ADD_ADVANCED_TEST( get_changed_trilinos_packages_test_${TEST_POSTFIX}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    OVERALL_NUM_MPI_PROCS 1
    ENVIRONMENT # We can't use real git so we use mockprogram.py instead!
      GCTP_GIT_OVERRIDE=${${PROJECT_NAME}_TRIBITS_DIR}/python_utils/mockprogram.py
      MOCKPROGRAM_INOUT_FILE_OVERRIDE=${ABS_TEST_DIR}/mockprogram_inout.txt

    TEST_0
      MESSAGE "Copy the mockprogram.py inout file each time since it gets modifed!"
      CMND "${CMAKE_COMMAND}"
        ARGS -E copy "${GIT_MOCKPROGRAM_TXT_FILE}" mockprogram_inout.txt

    TEST_1
      MESSAGE "Run get-changed-trilinos-packages.sh with mock git"
      CMND "${CMAKE_CURRENT_SOURCE_DIR}/get-changed-trilinos-packages.sh"
        ARGS sha1-from sha1-to packagesToEnable.cmake package_subproject_list.cmake
      PASS_REGULAR_EXPRESSION_ALL
        "CHANGED_PACKAGES_FULL_LIST='${CHANGED_PACKAGES_FULL_LIST}'"
	"CHANGED_PACKAGES_ST_LIST='${CHANGED_PACKAGES_ST_LIST}'"

    TEST_2 CMND ${CMAKE_COMMAND}
      ARGS -P packagesToEnable.cmake
      PASS_REGULAR_EXPRESSION_ALL "${EXPECTED_ENABLES_REGEX}"

    TEST_3 CMND cat
      ARGS package_subproject_list.cmake
      PASS_REGULAR_EXPRESSION_ALL
        "set[(]CTEST_LABELS_FOR_SUBPROJECTS ${EXPECTED_CTEST_LABELS_FOR_SUBPROJECTS_REGEX}[)]"

    )

ENDFUNCTION()


set(ALL_TOPLEVEL_PACKAGES_REGEX
  "TrilinosFrameworkTests TrilinosATDMConfigTests .*Kokkos KokkosKernels Teuchos .*Tpetra .*Thyra .*SEACAS .*Teko .*MueLu .*ShyLU_DD ShyLU .*Tempus Stokhos ROL .*Panzer Adelus TrilinosCouplings TrilinosBuildStats TrilinosInstallTests")

create_get_changed_trilinos_packages_test(ProjectsList_TeuchosCore
  "ProjectsList.cmake\npackages/teuchos/core/CMakeLists.txt"
  "ALL_PACKAGES,TeuchosCore"
  "ALL_PACKAGES,TeuchosCore"
  "Setting Trilinos_ENABLE_ALL_PACKAGES = ON;Setting Trilinos_ENABLE_TeuchosCore = ON"
  "${ALL_TOPLEVEL_PACKAGES_REGEX}"
  )
  # The above test ensures that it can detect global build files as well as
  # package files.  Above, we check several of the parent packages in Trilinos
  # and we want to ensure that no subpackages are listed in the array
  # CTEST_LABELS_FOR_SUBPROJECTS.  We also check to confirm that
  # PyTrilinos is excluded.  To do that, we list non EX
  # packages before and after them as they currently appear in the
  # Trilinos/PackagesList.cmake file.  If the PackagesList.cmake file changes,
  # then this test may have to change to accomidate that.  Therefore, we don't
  # want to list too many packages above in the regex.


create_get_changed_trilinos_packages_test(cmake_std_atdm_TeuchosCore_PyTrilinos_Compadre
  "cmake/std/atdm/anything\npackages/teuchos/core/CMakeLists.txt\npackages/PyTrilinos/anything\npackages/compadre/CMakeLists.txt"
  "TrilinosATDMConfigTests,TeuchosCore,Compadre"
  "TrilinosATDMConfigTests,TeuchosCore,Compadre"
  "Setting Trilinos_ENABLE_TrilinosATDMConfigTests = ON;Setting Trilinos_ENABLE_TeuchosCore = ON;Setting Trilinos_ENABLE_Compadre = ON"
  "${ALL_TOPLEVEL_PACKAGES_REGEX}"
  )
  # The above test ensures that cmake/std/atdm/ changes don't trigger global
  # builds and it makes sure that PT (TeuchosCore) and ST packages (Compadre) do
  # get enabled (but excluded ST package PyTrilinos does not get included).


create_get_changed_trilinos_packages_test(TeuchosCore_PyTrilinos2_Compadre
  "packages/teuchos/core/CMakeLists.txt\npackages/PyTrilinos2/anything\npackages/compadre/a.txt"
  "TeuchosCore,PyTrilinos2,Compadre"
  "TeuchosCore,Compadre"
  "Setting Trilinos_ENABLE_TeuchosCore = ON;Setting Trilinos_ENABLE_Compadre = ON"
  "${ALL_TOPLEVEL_PACKAGES_REGEX}"
  )
  # The above test ensures PT (TeuchosCore) and ST packages (Compadre) do get
  # enabled but EX packages (PyTrilinos2) don't.


create_get_changed_trilinos_packages_test(cmake_std_atdm
  "cmake/std/atdm/anything\ncmake/std/atdm/shiller/tweaks/INTEL-RELEASE-OPENMP-HSW.cmake"
  "TrilinosATDMConfigTests"
  "TrilinosATDMConfigTests"
  "Setting Trilinos_ENABLE_TrilinosATDMConfigTests = ON"
  "${ALL_TOPLEVEL_PACKAGES_REGEX}"
  )
  # The above test ensures changes to just cmake/std/atdm/ triggers the enable
  # of only the TrilinosATDMConfigTests package and does not trigger the
  # enable of every package.


create_get_changed_trilinos_packages_test(ShyLU_NodeHTS
  "packages/shylu/shylu_node/hts/test/CMakeLists.txt"
  "ShyLU_NodeHTS"
  "ShyLU_NodeHTS"
  "Setting Trilinos_ENABLE_ShyLU_NodeHTS = ON"
  "${ALL_TOPLEVEL_PACKAGES_REGEX}"
  )
  # The above test ensures changes to packages/shylu/shylu_node/hts/ trigger
  # testing of the ST ShyLU_NodeST subpackage.  (This will cause PRs like
  # #7229 to actually test code.)


# NOTE: Above, currently every top-level package is listed as a CDash
# subproject.  Once we can implement TriBITS forward-package logic, then the
# above tests will be checking a subset of packages.



TRIBITS_PACKAGE_POSTPROCESS()
