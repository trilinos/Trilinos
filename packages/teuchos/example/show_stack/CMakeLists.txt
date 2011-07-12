INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)

PACKAGE_ADD_EXECUTABLE(
  show_stack
  SOURCES cxx_main.cpp
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_STACKTRACE)
IF (${PACKAGE_NAME}_ENABLE_STACKTRACE)

  SET(PASS_REGEXS "")

  IF (${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS OR CMAKE_BUILD_TYPE STREQUAL DEBUG)
    APPEND_SET(PASS_REGEXS
      "Teuchos::show_stacktrace();"
      )
  ENDIF()

  IF (CMAKE_BUILD_TYPE STREQUAL DEBUG)
    APPEND_SET(PASS_REGEXS
      "teuchos/example/show_stack/cxx_main.cpp., line 66, in main()"
      "f()"
      "teuchos/example/show_stack/cxx_main.cpp., line 57, in f()"
      "g()"
      "teuchos/example/show_stack/cxx_main.cpp., line 50, in g()"
      "char .p = NULL. .p = 0"
      )
  ENDIF()

  
  PACKAGE_ADD_ADVANCED_TEST(
    show_stack
    TEST_0 EXEC show_stack
      PASS_REGULAR_EXPRESSION_ALL
        "Traceback .most recent call last.:"
        "Segfault caught. Printing stacktrace:"
        "Done. Exiting the program."
        ${PASS_REGEXS}
    OVERALL_NUM_MPI_PROCS 1
    )
  
    # NOTE: Above, I had to put in '.' for characters where regular expression
    # characters like '(', ')', and '*' would be miss interpreted.  However, I
    # think that this test is still very strong.

ENDIF()
