IF(${PACKAGE_NAME}_ENABLE_Belos)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)
  TRIBITS_ADD_EXECUTABLE(
    mf-example-01
    SOURCES mf_example_01.cpp
    COMM serial mpi
  )
ENDIF()
