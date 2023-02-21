
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/../disc-fe/test/equation_set)
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/../disc-fe/test/closure_model)
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/bcstrategy)

SET(main_driver_SOURCES
  main_driver.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  me_main_driver 
  SOURCES ${main_driver_SOURCES}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(me_main_driver_files
  SOURCE_FILES
    energy-ss.xml
    EXEDEPS me_main_driver
  )

TRIBITS_ADD_ADVANCED_TEST(
  me_main_driver_energy-ss 
  TEST_0 EXEC me_main_driver 
    ARGS  --i=energy-ss.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )
