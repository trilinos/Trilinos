INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#     amesos2
#     SOURCES amesos2.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu xpetra
#   )


# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#     simple1D-mergedSmoo
#     SOURCES mainMergedSmoo.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu xpetra
#   )

IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    simple1D-Emin
    SOURCES Emin.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu xpetra
    )
ENDIF ()
