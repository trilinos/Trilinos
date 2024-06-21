#TRIBITS_ADD_EXECUTABLE_AND_TEST(driver
#  SOURCES ${HEADERS} driver.cpp
#  NUM_MPI_PROCS 4
#  COMM mpi serial
#  ARGS "basic.xml"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_COPY_FILES_TO_BINARY_DIR(Zoltan2_copy_test_options
#  SOURCE_FILES basic.xml
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#  NOEXEPREFIX
#)

ADD_SUBDIRECTORY(driverinputs)

TRIBITS_ADD_EXECUTABLE(
   test_driver
   SOURCES ${HEADERS} test_driver.cpp
   COMM serial mpi
)


IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  TRIBITS_ADD_TEST(
    test_driver
    NAME parmetisSubComm
    NUM_MPI_PROCS 2
    COMM serial mpi
    ARGS
    "./driverinputs/parmetisSubComm.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

TRIBITS_ADD_TEST(
    test_driver
    NAME geometricTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/geometricTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME test578
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/test578.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  TRIBITS_ADD_TEST(
    test_driver
    NAME test7944
    COMM serial mpi
    ARGS
    "./driverinputs/test7944.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()


IF (${PACKAGE_NAME}_ENABLE_PuLP)
  TRIBITS_ADD_TEST(
    test_driver
    NAME pulpTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/pulpTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    test_driver
    NAME pulpTestNoDistrib
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/pulpTestNoDistrib.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    test_driver
    NAME pulpVwgtTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/pulpVwgtTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    test_driver
    NAME pulpVwgtTestNoDistrib
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/pulpVwgtTestNoDistrib.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    test_driver
    NAME xtrapulpVwgtTest2
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/xtrapulpVwgtTest2.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    test_driver
    NAME xtrapulpVwgtTest2NoDistrib
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/xtrapulpVwgtTest2NoDistrib.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

TRIBITS_ADD_TEST(
    test_driver
    NAME geometricVwgtTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/geometricVwgtTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME geometricVwgt2Test
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/geometricVwgt2Test.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

# MultiJagged Tests
TRIBITS_ADD_TEST(
    test_driver
    NAME multijaggedVwgt2Test
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/multijaggedVwgt2Test.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

IF (${PROJECT_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_TEST(
    test_driver
    NAME coordinateDriverTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/galeriTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

TRIBITS_ADD_TEST(
    test_driver
    NAME rcbDriverTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/rcbTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME mjDriverTest
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "./driverinputs/multiJaggedTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME mjTestPremigrate
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/multiJaggedPremigrateTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME zoltanWithGraphAdapter
    COMM serial mpi
    ARGS
    "./driverinputs/zoltanGraphAdapter.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

# Hypergraph Tests
IF(Zoltan2_ENABLE_HyperGraph)
TRIBITS_ADD_TEST(
    test_driver
    NAME zoltanWithGraphAdapterVwgt
    COMM serial mpi
    ARGS
    "./driverinputs/zoltanGraphAdapterVwgt.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME zoltanWithGraphAdapterEwgt
    COMM serial mpi
    ARGS
    "./driverinputs/zoltanGraphAdapterEwgt.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
    )
ENDIF()

TRIBITS_ADD_TEST(
    test_driver
    NAME chacoGraphMetricsTest
    COMM serial mpi
    ARGS
    "./driverinputs/chacoGraphMetricsTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME chacoGraphMetricsTestNoDistribute
    COMM serial mpi
    ARGS
    "./driverinputs/chacoGraphMetricsTestNoDistribute.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME chacoGraphMetricsTestNoGraph
    COMM serial mpi
    ARGS
    "./driverinputs/chacoGraphMetricsTestNoGraph.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME chacoGraphMetricsTestNoGraphNoDistribute
    COMM serial mpi
    ARGS
    "./driverinputs/chacoGraphMetricsTestNoGraphNoDistribute.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

IF (${PROJECT_NAME}_ENABLE_Pamgen)
  TRIBITS_ADD_TEST(
    test_driver
    NAME simplePamgenTest
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "./driverinputs/simplePamgen.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()


# Command that copies files to the executable directory.
# KDD: I am using this for Zoltan data files.  When Zoltan tests are not built,
# KDD: we can't point driver to files in build directory (as they won't be 
# KDD: there.  Also, we don't want to use hard-code paths in XML files.
# KDD: Copying the files here seems to be best, unless we can make CMAKE
# KDD: do macro substitution in the XML files.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_zoltan_simple_for_driver_tests_
  SOURCE_FILES
    simple.coords
    simple.graph
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../zoltan/test/ch_simple
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_zoltan_vwgt_for_driver_tests_
  SOURCE_FILES
    vwgt.coords
    vwgt.graph
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../zoltan/test/ch_vwgt
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_zoltan_vwgt2_for_driver_tests_
  SOURCE_FILES
    vwgt2.coords
    vwgt2.graph
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../zoltan/test/ch_vwgt2
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_zoltan_ewgt_for_driver_tests_
  SOURCE_FILES
    ewgt.coords
    ewgt.graph
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../zoltan/test/ch_ewgt
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_zoltan_nograph_for_driver_tests_
  SOURCE_FILES
    nograph.coords
    nograph.graph
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../zoltan/test/ch_nograph
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_zoltan_7944_for_driver_tests_
  SOURCE_FILES
    7944.graph
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../zoltan/test/ch_7944
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

