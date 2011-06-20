INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(main_driver_SOURCES
  main_driver.cpp
  user_app_NOXObserverFactory_Epetra.hpp
  user_app_NOXObserver_Epetra.hpp
  user_app_RythmosObserverFactory_Epetra.hpp
  user_app_RythmosObserver_Epetra.hpp
  Panzer_STK_ModelEvaluatorFactory_Epetra.hpp
  Panzer_STK_ModelEvaluatorFactory_EpetraT.hpp
  Panzer_STK_NOXObserverFactory_Epetra.hpp
  Panzer_STK_RythmosObserverFactory_Epetra.hpp
  )

PACKAGE_ADD_EXECUTABLE(
  main_driver
  SOURCES ${main_driver_SOURCES}
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(main_driver_files
  SOURCE_FILES energy-ss.xml energy-transient.xml
  EXEDEPS main_driver
  )

PACKAGE_ADD_ADVANCED_TEST(
  main_driver_energy-ss
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

PACKAGE_ADD_ADVANCED_TEST(
  main_driver_energy-transient
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-transient.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )