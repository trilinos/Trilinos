
# Required modules
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  exit_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp some_c_func.c 
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_C_EXCEPTIONS)
IF (${PACKAGE_NAME}_ENABLE_C_EXCEPTIONS)
  SET(PASS_CRITERIA STANDARD_PASS_OUTPUT)
ELSE()
  SET(PASS_CRITERIA PASS_REGULAR_EXPRESSION "some_c_func.c:[0-9]+: error code = 5")
ENDIF()

PACKAGE_ADD_TEST(
  exit_test  
  NUM_MPI_PROCS 1
  COMM serial mpi
  ${PASS_CRITERIA}
  )
