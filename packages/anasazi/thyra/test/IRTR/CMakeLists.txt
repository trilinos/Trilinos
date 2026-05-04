
IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_IRTR_test
    SOURCES cxx_main.cpp
    ARGS 
      "--verbose --hefty --sort=SR"
      "--verbose --skinny --sort=SR"
      "--verbose --hefty --sort=LR"
      "--verbose --skinny --sort=LR" 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF()
