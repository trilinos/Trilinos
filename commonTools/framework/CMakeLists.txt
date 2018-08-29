TRIBITS_PACKAGE(TrilinosFrameworkTests)


TRIBITS_ADD_ADVANCED_TEST( ProjectCiFileChangeLogic_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND ${PYTHON_EXECUTABLE} 
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ProjectCiFileChangeLogic_UnitTests.py -v
    PASS_REGULAR_EXPRESSION "OK"
  )


FUNCTION(create_get_changed_trilinos_packages_test  TEST_POSTFIX
  FILES_CHANGED
  CHANGED_PACKAGES_FULL_LIST
  CHANGED_PACKAGES_PT_LIST
  EXPECTED_ENABLES_REGEX
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
        ARGS sha1-from sha1-to packagesToEnable.cmake
      PASS_REGULAR_EXPRESSION_ALL
        "CHANGED_PACKAGES_FULL_LIST='${CHANGED_PACKAGES_FULL_LIST}'"
	"CHANGED_PACKAGES_PT_LIST='${CHANGED_PACKAGES_PT_LIST}'"

    TEST_2 CMND ${CMAKE_COMMAND}
      ARGS -P packagesToEnable.cmake
      PASS_REGULAR_EXPRESSION_ALL "${EXPECTED_ENABLES_REGEX}"

    )

ENDFUNCTION()


create_get_changed_trilinos_packages_test(ProjectsList_TeuchosCore
  "ProjectsList.cmake\npackages/teuchos/core/CMakeLists.txt"
  "ALL_PACKAGES,TeuchosCore"
  "ALL_PACKAGES,TeuchosCore"
  "Setting Trilinos_ENABLE_ALL_PACKAGES = ON;Setting Trilinos_ENABLE_TeuchosCore = ON"
  )
  # The above test ensures that it can detect global build files as well as
  # package files.


create_get_changed_trilinos_packages_test(cmake_std_atdm_TeuchosCore_PyTrilinos
  "cmake/std/atdm/anything\npackages/teuchos/core/CMakeLists.txt\npackages/PyTrilinos/anything"
  "TeuchosCore,PyTrilinos"
  "TeuchosCore"
  "Setting Trilinos_ENABLE_TeuchosCore = ON"
  )
  # The above test ensures that cmake/std/atdm/ changes don't trigger global
  # builds and it makes sure that ST packages (PyTrilinos) don't get enabled.


create_get_changed_trilinos_packages_test(cmake_std_atdm
  "cmake/std/atdm/anything\ncmake/std/atdm/shiller/tweaks/INTEL-RELEASE-OPENMP-HSW.cmake"
  ""
  ""
  "-- NOTE: No packages are being enabled!"
  )
  # The above test ensures that changes that don't anything are handled
  # smoothly


TRIBITS_PACKAGE_POSTPROCESS()
