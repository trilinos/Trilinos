INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/model_factory)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/bcstrategy)

PACKAGE_ADD_EXECUTABLE(
  main_driver
  SOURCES Panzer_MainDriver.cpp ${UNIT_TEST_DRIVER}
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(main_driver_files
  SOURCE_FILES energy-ss.xml energy-transient.xml
  EXEDEPS main_driver
  )

PACKAGE_ADD_ADVANCED_TEST(
  main_driver_energy
  TEST_0 EXEC main_driver 
    ARGS  --i=energy-ss.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  TEST_1 EXEC main_driver 
    ARGS  --i=energy-transient.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )