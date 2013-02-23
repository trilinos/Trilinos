INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

IF(NOT HAVE_MUELU_COMPLEX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Permutation
    SOURCES Permutation.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu xpetra
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Permutation_cp
SOURCE_FILES test1.txt test2.txt test3.txt test4.txt test5.txt exp1.txt exp2.txt exp3.txt exp4.txt exp5.txt
)
