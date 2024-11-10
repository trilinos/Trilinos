
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/closure_model)
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/bcstrategy)

SET(main_driver_SOURCES
  main_driver.cpp
  user_app_NOXObserverFactory.hpp
  user_app_NOXObserver_WriteToExodus.hpp
  user_app_NOXObserver_NeumannBCAnalyticSystemTest.hpp
  )

TRIBITS_ADD_EXECUTABLE(
  main_driver
  SOURCES ${main_driver_SOURCES}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(main_driver_files
  SOURCE_FILES
    energy-ss.xml
    energy-ss-tp.xml
    energy-ss-point-calc.xml
    energy-ss-tp-delay-prec.xml
    energy-ss-loca-eigenvalue.xml
    energy-ss-blocked.xml
    energy-transient-tempus-blocked.xml
    energy-neumann.xml
    energy-ss-blocked-tp.xml
    periodic_wedge.xml
    periodic_wedge.pam
    energy-periodic-32bit-limit.xml
  EXEDEPS main_driver
  )

IF(PANZER_HAVE_EPETRA_STACK)
  TRIBITS_ADD_ADVANCED_TEST(
    main_driver_energy-ss
    TEST_0 EXEC main_driver
      ARGS  --i=energy-ss.xml --exodus-io-num-procs=1
      PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
      XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
    )

  TRIBITS_ADD_ADVANCED_TEST(
    main_driver_energy-ss-blocked
    TEST_0 EXEC main_driver
      ARGS  --i=energy-ss-blocked.xml
      PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    )
ENDIF(PANZER_HAVE_EPETRA_STACK)

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-tp
  TEST_0 EXEC main_driver
    ARGS  --i=energy-ss-tp.xml --exodus-io-num-procs=1
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

  TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-point-calc
  TEST_0 EXEC main_driver
    ARGS  --i=energy-ss-point-calc.xml --exodus-io-num-procs=1 --point-calc
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    NUM_MPI_PROCS 4
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-tp-delay-prec
  TEST_0 EXEC main_driver
    ARGS  --i=energy-ss-tp-delay-prec.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-loca-eigenvalue
  TEST_0 EXEC main_driver
    ARGS  --i=energy-ss-loca-eigenvalue.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-ss-blocked-tp
  TEST_0 EXEC main_driver
    ARGS  --i=energy-ss-blocked-tp.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-periodic-32bit-limit
  TEST_0 EXEC main_driver
    ARGS  --i=energy-periodic-32bit-limit.xml
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
  )

TRIBITS_ADD_ADVANCED_TEST(
  main_driver_energy-neumann
  TEST_0 EXEC main_driver
    ARGS  --i=energy-neumann.xml --flux-calc
    PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    NUM_MPI_PROCS 1
  )

IF(${PACKAGE_NAME}_ENABLE_Teko)
  IF(${PACKAGE_NAME}_ENABLE_Tempus)
    TRIBITS_ADD_ADVANCED_TEST(
      main_driver_energy-transient-tempus-blocked
      TEST_0 EXEC main_driver
        ARGS  --i=energy-transient-tempus-blocked.xml
        PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
      )
  ENDIF()
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Pamgen)
  TRIBITS_ADD_ADVANCED_TEST(
    periodic_wedge
    TEST_0 EXEC main_driver
      ARGS  --i=periodic_wedge.xml
      PASS_REGULAR_EXPRESSION "panzer::MainDriver run completed."
    )
ENDIF()
