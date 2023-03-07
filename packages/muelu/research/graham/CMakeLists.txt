IF(${PACKAGE_NAME}_ENABLE_Belos)
  TRIBITS_ADD_EXECUTABLE(
    mf-example-01
    SOURCES mf_example_01.cpp
    COMM serial mpi
    )
ENDIF()
