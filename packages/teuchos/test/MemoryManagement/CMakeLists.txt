
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageAddAdvancedTest)


#
# The following tests are expressed in this order in order to catch compunding
# errors.  That is, if the ealier test executables fail, the other test
# executables are also likely to fail as well.
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Boost)
IF (${PACKAGE_NAME}_ENABLE_Boost)
  ASSERT_DEFINED(HAS_TEUCHOS_BOOST_IS_POLYMORPHIC)
  IF (HAS_TEUCHOS_BOOST_IS_POLYMORPHIC)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      GetBaseObjVoidPtr_UnitTests
      SOURCES GetBaseObjVoidPtr_UnitTests.cpp
        TestClasses.cpp
        ../UnitTest/Teuchos_StandardUnitTestMain.cpp
      NUM_MPI_PROCS 1
      ARGS --teuchos-suppress-startup-banner
      STANDARD_PASS_OUTPUT
      )
  ENDIF()
ENDIF()


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TestForException_throwTest
  SOURCES
    TestForException_throwTest.cpp
    SimpleThrowFunctions.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "This is an exception I throw"
  )


PACKAGE_ADD_ADVANCED_TEST(
  TestForException_throwTest_nostacktrace
  TEST_0 EXEC TestForException_throwTest
    ARGS  --no-enable-stacktrace
    PASS_REGULAR_EXPRESSION_ALL
      "teuchos/test/MemoryManagement/SimpleThrowFunctions.cpp:50"
      "Throw number = 1"
      "Throw test that evaluated to true: true"
      "This is an exception I throw!"
  OVERALL_NUM_MPI_PROCS 1
  )


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_STACKTRACE)
IF (${PACKAGE_NAME}_ENABLE_STACKTRACE)

  SET(PASS_REGEXS
    "teuchos/test/MemoryManagement/SimpleThrowFunctions.cpp:50"
    "Throw number = 1"
    "Throw test that evaluated to true: true"
    "This is an exception I throw!"
    )

  IF (CMAKE_BUILD_TYPE STREQUAL DEBUG)
    APPEND_SET(PASS_REGEXS
      "teuchos/test/MemoryManagement/TestForException_throwTest.cpp., line .+, in func_b()"
      "teuchos/test/MemoryManagement/TestForException_throwTest.cpp., line .+, in func_a()"
      )
  ELSEIF (${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS)
    APPEND_SET(PASS_REGEXS
      "teuchos/test/MemoryManagement/TestForException_throwTest.cpp., line .+, in func_.()"
      )
  ENDIF()

  PACKAGE_ADD_ADVANCED_TEST(
    TestForException_throwTest_stacktrace
    TEST_0 EXEC TestForException_throwTest
      ARGS  --enable-stacktrace
      PASS_REGULAR_EXPRESSION_ALL ${PASS_REGEXS}
    OVERALL_NUM_MPI_PROCS 1
    )
    # ToDo: Add more greps to look for the stacktrace output!

ENDIF()


# These RCPNodeHandle tests are pulled out to run first in order to catch
# problems before running the other memory management tests.
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RCPNodeHandle_UnitTests
  SOURCES
    RCPNodeHandle_UnitTests.cpp
    TestClasses.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MemoryManagement_UnitTests
  SOURCES
    RCP_UnitTests.cpp
    RCP_ForwardDeclUnitTests.cpp
    Ptr_UnitTests.cpp
    Array_UnitTest_helpers.cpp
    Array_UnitTests.cpp
    ArrayView_UnitTests.cpp
    ArrayRCP_UnitTests.cpp
    ArrayConversions_UnitTests.cpp
    TestClasses.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_DEBUG)

IF (${PACKAGE_NAME}_ENABLE_Boost AND HAS_TEUCHOS_BOOST_IS_POLYMORPHIC
  AND ${PACKAGE_NAME}_ENABLE_DEBUG
  )

  # These tests just provide additional verification that the unit tests show
  # that the RCPNode tracing code that catches tricky errors that the old
  # system could not catch.

  PACKAGE_ADD_TEST(
      MemoryManagement_UnitTests
      POSTFIX_AND_ARGS_0 duplicate_rcp_owning
        --test=duplicate_rcp_owning_polymorphic
        --details=ALL
      NUM_MPI_PROCS 1
      PASS_REGULAR_EXPRESSION
        "throws .*Teuchos.*DuplicateOwningRCPError.*: passed"
      )

  PACKAGE_ADD_TEST(
      MemoryManagement_UnitTests
      POSTFIX_AND_ARGS_0 duplicate_rcp_nonowning
        --test=rcp_duplicate_rcp_nonowning_polymorphic_different_addr
        --details=ALL
      NUM_MPI_PROCS 1
      PASS_REGULAR_EXPRESSION
        "throws .*Teuchos.*DanglingReferenceError.*: passed"
      )

  # NOTE: Above, we have to use the regular expressions in order to account
  # for machines that don't have name demangling turned on.  For example, on
  # GCC without name demangling it returns N7Teuchos22DanglingReferenceErrorE
  # instead of Teuchos::DanglingReferenceError.

ENDIF()


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleObjTbl_UnitTests
  SOURCES
    SimpleObjTbl_UnitTests.cpp
    TestClasses.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ArrayRCP_test  
  SOURCES ArrayRCP_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ArrayView_test  
  SOURCES ArrayView_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Array_test  
  SOURCES Array_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Ptr_test  
  SOURCES Ptr_test.cpp TestClasses.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RCP_test  
  SOURCES RCP_test.cpp TestClasses.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tuple_test  
  SOURCES Tuple_test.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  GCC_CheckedSTL_UnitTests
#  SOURCES GCC_CheckedSTL_UnitTests.cpp  
#    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS --teuchos-suppress-startup-banner
#  PASS_REGULAR_EXPRESSION
#    "error: attempt to dereference a past-the-end iterator"
#  )


#
# RCP Beginners Guide example code
#
IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  # The raw pointer version should build but is supposed to segfault so we will
  # not bother running it.  The problem is that even if we put WILL_FAIL onto
  # it, CTest will intercept signall 11 (segfault) and mark it as failed, even
  # though it is supposed to fail.  So much for WILL_FAIL ...
  PACKAGE_ADD_EXECUTABLE(
    example_raw_cpp_pointers
    SOURCES
      ${TRILINOS_HOME_DIR}/doc/RefCountPtr/example/example_raw_cpp_pointers.cpp
    )

  # The RCP version should run just fine
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    example_ref_count_ptr
    SOURCES
      ${TRILINOS_HOME_DIR}/doc/RefCountPtr/example/example_ref_count_ptr.cpp
    NUM_MPI_PROCS 1
    )
ENDIF()






#############################
# Performance tests
#############################


#
# RCP_PerforamnceTests
#


# Build the executable no matter what
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RCP_PerformanceTests
  NAME_POSTFIX basic
  CATEGORIES BASIC PERFORMANCE
  SOURCES
    RCP_Performance_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
  RCP_PerformanceTests
  CATEGORIES PERFORMANCE
  COMM serial mpi
  NUM_MPI_PROCS 1
  POSTFIX_AND_ARGS_0 performance
    --rel-cpu-speed=${${PROJECT_NAME}_REL_CPU_SPEED}
    --show-test-details=ALL --max-array-size=100000
    --max-rcp-create-destroy-ratio=1.05
    --max-rcp-raw-adjust-ref-count-ratio=20.0 --max-rcp-sp-adjust-ref-count-ratio=1.7
    --max-rcp-raw-obj-access-ratio=1.02
  STANDARD_PASS_OUTPUT
  )

#
# Array_PerforamnceTests
#

# Build the executable no matter what
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Array_PerformanceTests
  NAME_POSTFIX basic
  CATEGORIES BASIC PERFORMANCE
  SOURCES
    Array_Performance_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
  Array_PerformanceTests
  CATEGORIES PERFORMANCE
  COMM serial mpi
  NUM_MPI_PROCS 1
  POSTFIX_AND_ARGS_0 performance
    --rel-cpu-speed=${${PROJECT_NAME}_REL_CPU_SPEED}
    --details=ALL --max-array-size=100000
    --max-array-bracket-ratio=1.02 --max-array-iter-ratio=1.02
    --max-arrayrcp-self-iter-ratio=1.5
  STANDARD_PASS_OUTPUT
  )

