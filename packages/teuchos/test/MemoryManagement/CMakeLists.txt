
INCLUDE(PackageAddExecutableAndTest)


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
