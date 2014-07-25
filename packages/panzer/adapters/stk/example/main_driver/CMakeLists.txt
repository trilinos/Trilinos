INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddAdvancedTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(main_driver_SOURCES
  main_driver.cpp
  user_app_NOXObserverFactory.hpp
  user_app_NOXObserver_WriteToExodus.hpp
  user_app_NOXObserver_NeumannBCAnalyticSystemTest.hpp
  user_app_RythmosObserverFactory.hpp
  user_app_RythmosObserver_WriteToExodus.hpp
  )

TRIBITS_ADD_EXECUTABLE(
  main_driver
  SOURCES ${main_driver_SOURCES}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(main_driver_files
  SOURCE_FILES
    energy-ss.xml
    energy-ss-tp.xml
    energy-ss-loca-eigenvalue.xml
    energy-transient.xml
    energy-ss-blocked.xml
    energy-transient-blocked.xml
    energy-neumann.xml
    energy-ss-blocked-tp.xml
meshmotion.xml
  EXEDEPS main_driver
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss --exodus-io-num-procs=1
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  TEST_1 EXEC main_driver 
    ARGS  --i=energy-ss-tp.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-transient
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-transient.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-blocked
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss-blocked.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-loca-eigenvalue
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss-loca-eigenvalue.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-blocked-tp
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss-blocked-tp.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-neumann
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-neumann.xml --flux-calc
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    NUM_MPI_PROCS 1
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_meshmotion
  TEST_0 EXEC main_driver 
    ARGS  --i=meshmotion.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    NUM_MPI_PROCS 2
  )

IF(${PACKAGE_NAME}_ENABLE_Teko)
  TRIBITS_ADD_ADVANCED_TEST(
    main_driver_energy-transient-blocked
    TEST_0 EXEC main_driver 
      ARGS  --i=energy-transient-blocked.xml
      PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
      XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
    )
ENDIF()
