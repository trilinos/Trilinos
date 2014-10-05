INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Ifpack
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  )

IF ((${PACKAGE_NAME}_ENABLE_Amesos AND ${PACKAGE_NAME}_ENABLE_Ifpack) OR (${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Ifpack2))

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultipleFactoryManager
    SOURCES MultipleFactoryManager.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
    )

ENDIF()
