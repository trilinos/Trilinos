
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/zz)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/lb)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/Utilities/shared)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/params)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/order)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/par)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/tpls)

#
# Test of HIER_ASSIST option of hierarchical partitioning
#

IF (NOT ${PACKAGE_NAME}_ENABLE_PURIFY)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    stressTestColor
    SOURCES stressTestColor.c
    COMM mpi 
    NUM_MPI_PROCS 8
    FAIL_REGULAR_EXPRESSION "FAIL"
    ARGS
      "--size=10000"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS AND ${PACKAGE_NAME}_ENABLE_Scotch)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH_P_S
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    FAIL_REGULAR_EXPRESSION "FAIL"
    ARGS
      "--topology=2,4 --size=10000 --graph_package=parmetis"
      "--topology=2,4 --size=10000 --graph_package=scotch"
      "--topology=2,2,2 --size=10000 --graph_package=parmetis"
      "--topology=2,2,2 --size=10000 --graph_package=scotch"
  )

ELSEIF (${PACKAGE_NAME}_ENABLE_ParMETIS)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH_P
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    FAIL_REGULAR_EXPRESSION "FAIL"
    ARGS
      "--topology=2,4 --size=10000 --graph_package=parmetis"
      "--topology=2,2,2 --size=10000 --graph_package=parmetis"
  )

ELSEIF (${PACKAGE_NAME}_ENABLE_Scotch)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH_S
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    FAIL_REGULAR_EXPRESSION "FAIL"
    ARGS
      "--topology=2,4 --size=10000 --graph_package=scotch"
      "--topology=2,2,2 --size=10000 --graph_package=scotch"
  )

ENDIF()

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  stressTestGRAPH
  SOURCES stressTestGRAPH.c 
  COMM mpi 
  NUM_MPI_PROCS 8
  FAIL_REGULAR_EXPRESSION "FAIL"
  ARGS
    "--topology=2,4 --size=10000"
    "--topology=2,2,2 --size=10000"
)


#=============================================
#  Might be useful for scalbility tests
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  stressTestRCB
#  SOURCES stressTestRCB.c 
#  COMM mpi
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  ARGS
#    "10000"
#  )
#
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  stressTestRIB
#  SOURCES stressTestRIB.c 
#  COMM mpi
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  ARGS
#    "10000"
#  )
#=============================================

#
# Test of 64-bit global IDs (if built with 64 bit IDs)
#
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  stressTestPHG
  SOURCES stressTestPHG.c 
  COMM mpi 
  NUM_MPI_PROCS 8
  FAIL_REGULAR_EXPRESSION "FAIL"
  ARGS
    "--size=10000 --use_high_order_bits --use_varying_weights"
    "--use_graph_queries --size=10000 --use_high_order_bits --use_varying_weights"
  )
