INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(Zoltan/src/zz)

#
# Test of HIER_ASSIST option of hierarchical partitioning
#

IF (${PACKAGE_NAME}_ENABLE_ParMETIS AND ${PACKAGE_NAME}_ENABLE_Scotch)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH_P_S
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    ARGS
      "--topology=2,4 --size=10000 --graph_package=parmetis"
      "--topology=2,4 --size=10000 --graph_package=scotch"
      "--topology=2,2,2 --size=10000 --graph_package=parmetis"
      "--topology=2,2,2 --size=10000 --graph_package=scotch"
  )

ELSEIF (${PACKAGE_NAME}_ENABLE_ParMETIS)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH_P
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    ARGS
      "--topology=2,4 --size=10000 --graph_package=parmetis"
      "--topology=2,2,2 --size=10000 --graph_package=parmetis"
  )

ELSEIF (${PACKAGE_NAME}_ENABLE_Scotch)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH_S
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    ARGS
      "--topology=2,4 --size=10000 --graph_package=scotch"
      "--topology=2,2,2 --size=10000 --graph_package=scotch"
  )

ELSE()
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stressTestGRAPH
    SOURCES stressTestGRAPH.c 
    COMM mpi 
    NUM_MPI_PROCS 8
    ARGS
      "--topology=2,4 --size=10000"
      "--topology=2,2,2 --size=10000"
    )
ENDIF()


#=============================================
#  Might be useful for scalbility tests
#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  stressTestRCB
#  SOURCES stressTestRCB.c 
#  COMM mpi
#  ARGS
#    "10000"
#  )
#
#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  stressTestRIB
#  SOURCES stressTestRIB.c 
#  COMM mpi
#  ARGS
#    "10000"
#  )
#=============================================

#
# Test of 64-bit global IDs (if built with 64 bit IDs)
#
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  stressTestPHG
  SOURCES stressTestPHG.c 
  COMM mpi 
  NUM_MPI_PROCS 8
  ARGS
    "--size=10000 --use_high_order_bits --use_varying_weights"
    "--use_graph_queries --size=10000 --use_high_order_bits --use_varying_weights"
  )
