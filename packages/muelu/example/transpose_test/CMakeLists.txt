INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  )

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  transpose_test
  SOURCES transpose_test.cpp
  COMM mpi serial
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(transpose_test_example_cp
  SOURCE_FILES Test.mat
  )

ENDIF()
