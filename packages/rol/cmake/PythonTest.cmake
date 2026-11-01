MACRO(PYTHON_ADD_TEST TEST_NAME)

  TRIBITS_ADD_TEST(
    ${Python3_EXECUTABLE}
    NOEXEPREFIX
    NOEXESUFFIX
    NAME PYTHON_${TEST_NAME}
    ARGS "-B -m unittest ${TEST_NAME}.py"
    NUM_MPI_PROCS 1
    ${ARGN}
    )

  # The tests run from their source directory, so that they always see the
  # current files; -B keeps Python from writing bytecode there.
  IF(TPL_ENABLE_MPI)
    tribits_set_tests_properties(ROL_PYTHON_${TEST_NAME}_MPI_1
      PROPERTIES ENVIRONMENT "${PyROL_PYTHONPATH}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  ELSE()
    tribits_set_tests_properties(ROL_PYTHON_${TEST_NAME}
      PROPERTIES ENVIRONMENT "${PyROL_PYTHONPATH}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  ENDIF()

ENDMACRO(PYTHON_ADD_TEST TEST_NAME)
