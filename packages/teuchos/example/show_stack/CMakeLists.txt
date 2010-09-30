INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)

PACKAGE_ADD_EXECUTABLE(
  show_stack
  SOURCES cxx_main.cpp
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_STACKTRACE)
IF (${PACKAGE_NAME}_ENABLE_STACKTRACE)
  
  PACKAGE_ADD_ADVANCED_TEST(
    show_stack
    TEST_0 EXEC show_stack
      PASS_REGULAR_EXPRESSION_ALL
        "Traceback .most recent call last.:"
        "teuchos/example/show_stack/cxx_main.cpp., line 16, in main()"
        "f()"
        "teuchos/example/show_stack/cxx_main.cpp., line 10, in f()"
        "g()"
        "teuchos/example/show_stack/cxx_main.cpp., line 5, in g()"
        "Teuchos::show_stacktrace();"
        "teuchos/src/Teuchos_stacktrace.cpp., .+, in Teuchos::show_stacktrace()"
        "teuchos/src/Teuchos_stacktrace.cpp., .+, in Teuchos::get_stacktrace()"
        "Segfault caught. Printing stacktrace:"
        "teuchos/example/show_stack/cxx_main.cpp., line 19, in main()"
        "char .p = NULL. .p = 0"
    OVERALL_NUM_MPI_PROCS 1
    )
  
    # NOTE: Above, I had to put in '.' for characters where regular expression
    # characters like '(', ')', and '*' would be miss interpreted.  However, I
    # think that this test is still very strong.

ENDIF()
