INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddAdvancedTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(main_driver_SOURCES
  main_driver.cpp
  user_app_NOXObserverFactory_Epetra.hpp
  user_app_NOXObserver_EpetraToExodus.hpp
  user_app_NOXObserver_NeumannBCAnalyticSystemTest.hpp
  user_app_RythmosObserverFactory_Epetra.hpp
  user_app_RythmosObserver_EpetraToExodus.hpp
  Panzer_STK_ModelEvaluatorFactory_Epetra.hpp
  Panzer_STK_ModelEvaluatorFactory_Epetra_impl.hpp
  Panzer_STK_NOXObserverFactory.hpp
  Panzer_STK_RythmosObserverFactory.hpp
  )

TRIBITS_ADD_EXECUTABLE(
  main_driver
  SOURCES ${main_driver_SOURCES}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(main_driver_files
  SOURCE_FILES energy-ss.xml energy-transient.xml energy-ss-blocked.xml
energy-transient-blocked.xml energy-neumann.xml
  EXEDEPS main_driver
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss --exodus-io-num-procs=1
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-transient
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-transient.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-blocked
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss-blocked.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-neumann
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-neumann.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    NUM_MPI_PROCS 1
  )

IF(${PACKAGE_NAME}_ENABLE_Teko)
  TRIBITS_ADD_ADVANCED_TEST(
    main_driver_energy-transient-blocked
    TEST_0 EXEC main_driver 
      ARGS  --i=energy-transient-blocked.xml
      PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    )
ENDIF()
