TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

# GH: remove all the permutation tests based on this?
IF(NOT HAVE_MUELU_COMPLEX
   AND ${PACKAGE_NAME}_ENABLE_EpetraExt
   AND ${PACKAGE_NAME}_ENABLE_Epetra
   AND NOT ${PACKAGE_NAME}_ENABLE_Tpetra
)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Permutation
    SOURCES Permutation.cpp
    ARGS "--linAlgebra=Epetra"
    COMM serial mpi
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Permutation_cp
  SOURCE_FILES
    exp1.txt
    exp2.txt
    exp3.txt
    exp4.txt
    exp5.txt
    test1.txt
    test2.txt
    test3.txt
    test4.txt
    test5.txt
)
