
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    rcbTest
    SOURCES rcbTest.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()
  
# MultiJagged Tests
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    mjTest
    SOURCES MultiJaggedTest.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "I=1.03 P=2,2 O=0 C=4 F=simple TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_simple_rectilinear_boxes
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "I=1.03 P=2,2 O=0 C=4 F=simple R=1 TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "P=2,2,2,2 C=16 O=1 F=GeomGenParam.txt TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

#  TRIBITS_ADD_TEST(
#    mjTest
#    NAME mjTest_geomgen_boxes
#    NUM_MPI_PROCS 3
#    COMM serial mpi
#    ARGS
#    "P=2,2,2,2 C=16 O=1 F=GeomGenParam.txt TB=1"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#  )

TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_noholes
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "P=2,2,2,2 C=16 O=1 F=GeomGenParam_noholes.txt TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

#muelu like usecase, don't specify mjparts
TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_nomjparts
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "C=16 O=1 F=GeomGenParam_noholes.txt TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

#  TRIBITS_ADD_TEST(
#    mjTest
#    NAME mjTest_geomgen_nomjparts_boxes
#    NUM_MPI_PROCS 3
#    COMM serial mpi
#    ARGS
#    "C=16 O=1 F=GeomGenParam_noholes.txt TB=1"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#  )

#muelu like usecase, don't specify mjparts, odd parts
TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_nomjparts_odd
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "C=9 O=1 F=GeomGenParam_noholes.txt TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

#  TRIBITS_ADD_TEST(
#    mjTest
#    NAME mjTest_geomgen_nomjparts_odd_boxes
#    NUM_MPI_PROCS 3
#    COMM serial mpi
#    ARGS
#    "C=9 O=1 F=GeomGenParam_noholes.txt TB=1"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  taskMapper3
  SOURCES taskMappingTest3.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS
  "PROC=machine.txt NX=8 NY=8 NZ=8"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  blockTest
  SOURCES blockTest.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  pamgenMeshAdapterTest
  NAME pamgenMeshAdapterTest_color
  SOURCES pamgenMeshAdapterTest.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
    "--action=color"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  pamgenMeshAdapterTest
  NAME pamgenMeshAdapterTest_mj
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
    "--action=mj --nparts=6"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

IF (${PACKAGE_NAME}_ENABLE_Scotch)
  TRIBITS_ADD_TEST(
    pamgenMeshAdapterTest
    NAME pamgenMeshAdapterTest_scotch
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
      "--action=scotch --nparts=6"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Partitioning1
  SOURCES partitioning1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_OneProc
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_VWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --vertexWeights=3"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_OneProc_VWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute --vertexWeights=1"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_EWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --edgeWeights=1"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_OneProc_EWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute --edgeWeights=2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_ParMETIS
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_ParMETIS_OneProc
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute --method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_ParMETIS_VWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --vertexWeights=3 --method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_ParMETIS_OneProc_VWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute --vertexWeights=1 --method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_ParMETIS_EWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --edgeWeights=1 --method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Partitioning1
  NAME Partitioning1_ParMETIS_OneProc_EWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute --edgeWeights=2 --method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

##########################################
#IF (${PACKAGE_NAME}_ENABLE_Experimental AND ${PACKAGE_NAME}_ENABLE_Experimental_Wolf)

IF (${PACKAGE_NAME}_ENABLE_Experimental)

IF (${PACKAGE_NAME}_ENABLE_Experimental_Wolf)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Wolf
  SOURCES wolf.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()
#
ENDIF()
##########################################


# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_partition_tests
  SOURCE_FILES
    simple.mtx
    simple_coord.mtx
    USAir97.mtx
    USAir97_coord.mtx
    GeomGenParam.txt
    GeomGenParam_noholes.txt
    proc.txt
    part.txt
    machine.txt
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_more_files_for_partition_tests
  SOURCE_FILES
    Poisson.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

