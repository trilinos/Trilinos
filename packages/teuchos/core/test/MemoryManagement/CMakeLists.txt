#
# The following tests are expressed in this order in order to catch compunding
# errors.  That is, if the ealier test executables fail, the other test
# executables are also likely to fail as well.
#

# FIXME: Remove this hard-coded setting.  We always have C++11 now.
SET(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR  TRUE)

IF (HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    GetBaseObjVoidPtr_UnitTests
    SOURCES GetBaseObjVoidPtr_UnitTests.cpp
      TestClasses.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    NUM_MPI_PROCS 1
    ARGS --teuchos-suppress-startup-banner
    STANDARD_PASS_OUTPUT
    )
ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AssertAndThrow_UnitTests
  SOURCES
    AssertAndThrow_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TestForException_throwTest
  SOURCES
    TestForException_throwTest.cpp
    SimpleThrowFunctions.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "This is an exception I throw"
  )


TRIBITS_ADD_ADVANCED_TEST(
  TestForException_throwTest_nostacktrace
  TEST_0 EXEC TestForException_throwTest
    ARGS  --no-enable-stacktrace
    PASS_REGULAR_EXPRESSION_ALL
      "teuchos.core.test.MemoryManagement.SimpleThrowFunctions.cpp:16"
      "Throw number = 1"
      "Throw test that evaluated to true: true"
      "This is an exception I throw!"
  OVERALL_NUM_MPI_PROCS 1
  )


ASSERT_DEFINED(${PARENT_PACKAGE_NAME}_ENABLE_STACKTRACE)
IF (${PARENT_PACKAGE_NAME}_ENABLE_STACKTRACE)

  SET(PASS_REGEXS
    "teuchos.core.test.MemoryManagement.SimpleThrowFunctions.cpp:16"
    "Throw number = 1"
    "Throw test that evaluated to true: true"
    "This is an exception I throw!"
    )

  IF (CMAKE_BUILD_TYPE STREQUAL DEBUG)
    APPEND_SET(PASS_REGEXS
      "teuchos.core.test.MemoryManagement.TestForException_throwTest.cpp., line .+, in func_b()"
      "teuchos.core.test.MemoryManagement.TestForException_throwTest.cpp., line .+, in func_a()"
      )
  ELSEIF (${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS)
    APPEND_SET(PASS_REGEXS
      "teuchos.core.test.MemoryManagement.TestForException_throwTest.cpp., line .+, in func_.()"
      )
  ENDIF()

  TRIBITS_ADD_ADVANCED_TEST(
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
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RCPNodeHandle_UnitTests
  SOURCES
    RCPNodeHandle_UnitTests.cpp
    TestClasses.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RCPNodeTracer_UnitTests
  SOURCES
    RCPNodeTracer_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
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
    Any_UnitTests.cpp
    TestClasses.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MemoryManagement_MT_UnitTests
  EXCLUDE_IF_NOT_TRUE
    ${PARENT_PACKAGE_NAME}_ENABLE_THREAD_SAFE_UNIT_TESTS
  SOURCES
    Ptr_MT_UnitTests.cpp
    ArrayView_MT_UnitTests.cpp
    Array_MT_UnitTests.cpp
    Tuple_MT_UnitTests.cpp
    ArrayRCP_MT_UnitTests.cpp
    RCP_MT_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  NUM_TOTAL_CORES_USED 4
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MemoryManagement_MT_Nightly_UnitTests
  EXCLUDE_IF_NOT_TRUE
    ${PARENT_PACKAGE_NAME}_ENABLE_THREAD_SAFE_UNIT_TESTS
  SOURCES
    Ptr_MT_Nightly_UnitTests.cpp
    ArrayView_MT_Nightly_UnitTests.cpp
    Array_MT_Nightly_UnitTests.cpp
    Tuple_MT_Nightly_UnitTests.cpp
    ArrayRCP_MT_Nightly_UnitTests.cpp
    RCP_MT_Nightly_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  CATEGORIES NIGHTLY
  NUM_MPI_PROCS 1
  NUM_TOTAL_CORES_USED 4
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


ASSERT_DEFINED(${PARENT_PACKAGE_NAME}_ENABLE_DEBUG)

IF (${PARENT_PACKAGE_NAME}_ENABLE_DEBUG)
  TRIBITS_ADD_EXECUTABLE(
    MemoryManagement_RCP_Abort
    SOURCES
      RCP_Abort.cpp
    )

  TRIBITS_ADD_ADVANCED_TEST(
    MemoryManagement_RCP_Abort_Verify
    TEST_0 EXEC MemoryManagement_RCP_Abort
      ARGS  --enable-stacktrace
      PASS_REGULAR_EXPRESSION_ALL
        "PROGRAM ABORTING"
    OVERALL_NUM_MPI_PROCS 1
    )
ENDIF()


IF (HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR  AND  ${PARENT_PACKAGE_NAME}_ENABLE_DEBUG)

  # These tests just provide additional verification that the unit tests show
  # that the RCPNode tracing code that catches tricky errors that the old
  # system could not catch.

  TRIBITS_ADD_TEST(
      MemoryManagement_UnitTests
      POSTFIX_AND_ARGS_0 duplicate_rcp_owning
        --test=duplicate_rcp_owning_polymorphic
        --details=ALL
      NUM_MPI_PROCS 1
      PASS_REGULAR_EXPRESSION
        "throws .*Teuchos.*DuplicateOwningRCPError.*: passed"
      )

  TRIBITS_ADD_TEST(
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


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SimpleObjTbl_UnitTests
  SOURCES
    SimpleObjTbl_UnitTests.cpp
    TestClasses.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MemoryManagementExt_UnitTests
  SOURCES
    SimpleObjectDB_UnitTests.cpp
    TestClasses.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ArrayRCP_test  
  SOURCES ArrayRCP_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ArrayView_test  
  SOURCES ArrayView_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Array_test  
  SOURCES Array_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Ptr_test  
  SOURCES Ptr_test.cpp TestClasses.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RCP_test  
  SOURCES RCP_test.cpp TestClasses.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tuple_test  
  SOURCES Tuple_test.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  GCC_CheckedSTL_UnitTests
#  SOURCES GCC_CheckedSTL_UnitTests.cpp  
#    ${TEUCHOS_STD_UNIT_TEST_MAIN}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS --teuchos-suppress-startup-banner
#  PASS_REGULAR_EXPRESSION
#    "error: attempt to dereference a past-the-end iterator"
#  )


#
# Test Teuchos::RCPNodeTracer::setPrintActiveRcpNodesOnExit()
#
# Writing a test case for this feature is hard because the only way to test it
# is to run a full program.  Therefore, this is not really a unit test but it
# is still an effective verification test.
#

FUNCTION(TEUCHOS_TEST_setPrintActiveRcpNodesOnExit)

  # Need a function to create a private scope so we can turne these variables
  # off.
  SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF)
  SET(${PROJECT_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF)

  TRIBITS_ADD_EXECUTABLE(
    testDisablePrintActiveRcpNodesOnExit
    SOURCES testDisablePrintActiveRcpNodesOnExit.cpp TestClasses.cpp
    )

  TRIBITS_ADD_ADVANCED_TEST(
    testDisablePrintActiveRcpNodesOnExit_print
    TEST_0
      EXEC testDisablePrintActiveRcpNodesOnExit 
      NUM_MPI_PROCS 1
      PASS_REGULAR_EXPRESSION_ALL
        "End Result: TEST PASSED"
        "The following Teuchos::RCPNode objects were created"
    )
  # ABOVE: If the you accidentally disable the printing of the active RCP
  # nodes, the second regex will not match and the test will fail.  Try it for
  # yourself!

  TRIBITS_ADD_TEST(
    testDisablePrintActiveRcpNodesOnExit
    POSTFIX_AND_ARGS_0  noprint  --no-print-active-rcp-nodes-on-exit
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION
      "End Result: TEST PASSED"
    FAIL_REGULAR_EXPRESSION
      "The following Teuchos::RCPNode objects were created"
    )
  # ABOVE: If you the above flag --no-print-active-rcp-nodes-on-exit does not
  # set the bool to turn off the printing, then FAIL_REGULAR_EXRESSION clause
  # will kick in and cause the test to fail.  Try it for yourself.

ENDFUNCTION()

IF (${PARENT_PACKAGE_NAME}_ENABLE_DEBUG)
  TEUCHOS_TEST_setPrintActiveRcpNodesOnExit()
ENDIF()


#
# Test TEUCHOS_TEST_FOR_TERMINATION()
#x1
# This macro aborts the program so the only way to test it is to run a program.
#

TRIBITS_ADD_EXECUTABLE( testTeuchosTestForTermination
  SOURCES testTeuchosTestForTermination.cpp
  )
TRIBITS_ADD_ADVANCED_TEST( testTeuchosTestForTermination
  TEST_0
    EXEC testTeuchosTestForTermination
    ARGS --teuchos-suppress-startup-banner  # Avoid manging multi-proc output
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION
      "Bingo, we are terminating on procid == terminate_on_procid = 0"
    ALWAYS_FAIL_ON_ZERO_RETURN
  )
# NOTE: Above, we have to use TRIBITS_ADD_ADVANCED_TEST() and not
# TRIBITS_ADD_TEST() to to grep for the expected output which is made to
# std::cerr and not std::cout.  This test fails with a serial build if you use
# TRIBITS_ADD_TEST() which just uses raw ADD_TEST() and the CTest property
# PASS_REGULAR_EXPRESSION.  It appears that at least some versions of CTest
# will only grep grep STDOUT and not STDERR when applying
# PASS_REGULAR_EXPRESSION.  But with an MPI build, `mpirun` appears to combine
# STDOUT and STDERR into a sinlge STDOUT somehow that allows CTest to find the
# grep.  Using TRIBITS_ADD_ADVANCED_TEST() solves the problem because it uses
# EXECUTE_PROCESS() to tie STDOUT and STDERR together and then its own
# implementation of PASS_REGULAR_EXPRESSION finds the text send to STDERR.

SET(testTeuchosTestForTermination_NP 4)
MATH(EXPR testTeuchosTestForTermination_NP_m_1
  "${testTeuchosTestForTermination_NP}-1")
FOREACH(PID RANGE ${testTeuchosTestForTermination_NP_m_1})
  TRIBITS_ADD_ADVANCED_TEST(
    testTeuchosTestForTermination_${PID}_MPI_${testTeuchosTestForTermination_NP}
    COMM mpi
    TEST_0
      MESSAGE "Call TEUCHOS_TEST_FOR_TERMINATION(true) only on PID=${PID} which should terminate all processes"
      EXEC testTeuchosTestForTermination
      ARGS --terminate-on-procid=${PID}
        --teuchos-suppress-startup-banner  # Avoid manging multi-proc output
      NUM_MPI_PROCS ${testTeuchosTestForTermination_NP}
      PASS_REGULAR_EXPRESSION_ALL
        "p=${PID}: .*/testTeuchosTestForTermination.cpp"
        "Bingo, we are terminating on procid == terminate_on_procid = ${PID}"
      ALWAYS_FAIL_ON_ZERO_RETURN
    )
  # NOTE: If each of these tests does not terminate right away, then that is a
  # sign that std::terminate() on a single MPI process does not cause all other
  # processes to also terminate.
ENDFOREACH()


# Test call to terminate on bad destructor call
IF (${PARENT_PACKAGE_NAME}_ENABLE_DEBUG)
  TRIBITS_ADD_EXECUTABLE( testTeuchosRCPNodeImplDestructTerminate
    SOURCES testTeuchosRCPNodeImplDestructTerminate.cpp
    )
  TRIBITS_ADD_ADVANCED_TEST( testTeuchosRCPNodeImplDestructTerminate
    TEST_0
      EXEC testTeuchosRCPNodeImplDestructTerminate
      ARGS --teuchos-suppress-startup-banner  # Avoid manging multi-proc output
      NUM_MPI_PROCS 1
      PASS_REGULAR_EXPRESSION
        "Error, the underlying object must be explicitly deleted before deleting the node object"
      ALWAYS_FAIL_ON_ZERO_RETURN
    )
  # NOTE: Above, we have to use TRIBITS_ADD_ADVANCED_TEST() and not
  # TRIBITS_ADD_TEST() to grep the output with PASS_REGULAR_EXPRESSION for the
  # same reason as for the test testTeuchosTestForTermination above.
ENDIF()


#
# RCP Beginners Guide example code
#
SET(EXAMPLE_RAW_CPP_POINTERS_CPP
  ${Trilinos_SOURCE_DIR}/doc/RefCountPtr/example/example_raw_cpp_pointers.cpp)
IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE AND EXISTS ${EXAMPLE_RAW_CPP_POINTERS_CPP})
  # The raw pointer version should build but is supposed to segfault so we will
  # not bother running it.  The problem is that even if we put WILL_FAIL onto
  # it, CTest will intercept signall 11 (segfault) and mark it as failed, even
  # though it is supposed to fail.  So much for WILL_FAIL ...
  TRIBITS_ADD_EXECUTABLE(
    example_raw_cpp_pointers
    SOURCES
      ${Trilinos_SOURCE_DIR}/doc/RefCountPtr/example/example_raw_cpp_pointers.cpp
    )

  # The RCP version should run just fine
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_ref_count_ptr
    SOURCES
      ${Trilinos_SOURCE_DIR}/doc/RefCountPtr/example/example_ref_count_ptr.cpp
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
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RCP_PerformanceTests
  NAME_POSTFIX basic
  CATEGORIES BASIC PERFORMANCE
  SOURCES
    RCP_Performance_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
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
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Array_PerformanceTests
  NAME_POSTFIX basic
  CATEGORIES BASIC PERFORMANCE
  SOURCES
    Array_Performance_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
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

