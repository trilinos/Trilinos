
IF (${PROJECT_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
      rcbTest
      SOURCES rcbTest.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )
ENDIF()
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    fix2010
    SOURCES fix2010.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

# MultiJagged Tests
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    mj_int_coordinates
    SOURCES mj_int_coordinates.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    mj_backwardcompat
    SOURCES mj_backwardcompat.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    fix4785
    SOURCES fix4785.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

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

IF (Tpetra_ENABLE_CUDA)

TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_uvmOff
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "I=1.03 P=2,2 O=0 C=4 F=simple TB=1 UVM=0"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_uvmOff
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "P=2,2,2,2 C=16 O=1 F=GeomGenParam.txt TB=1 UVM=0"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_noholes_uvmOff
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "P=2,2,2,2 C=16 O=1 F=GeomGenParam_noholes.txt TB=1 UVM=0"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

#muelu like usecase, don't specify mjparts
TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_nomjparts_uvmOff
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "C=16 O=1 F=GeomGenParam_noholes.txt TB=1 UVM=0"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

#muelu like usecase, don't specify mjparts, odd parts
TRIBITS_ADD_TEST(
    mjTest
    NAME mjTest_geomgen_nomjparts_odd_uvmOff
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "C=9 O=1 F=GeomGenParam_noholes.txt TB=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TaskMappingTest3
  SOURCES TaskMappingTest3.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS
  "PROC=machine.txt NX=8 NY=8 NZ=8"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TaskMappingTest
  SOURCES TaskMappingTest.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "nx 8 ny 8 nz 8"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TaskMappingProblemTest
  SOURCES TaskMappingProblemTest.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "nx 8 ny 8 nz 8"
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

IF (${PROJECT_NAME}_ENABLE_Pamgen)

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

  TRIBITS_ADD_TEST(
    pamgenMeshAdapterTest
    NAME pamgenMeshAdapterTest_zoltanrcb
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
      "--action=zoltan_rcb --nparts=6"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    pamgenMeshAdapterTest
    NAME pamgenMeshAdapterTest_zoltanhg
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
      "--action=zoltan_hg --nparts=6"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    pamgenMeshAdapterTest
    NAME pamgenMeshAdapterTest_ghosthg
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
      "--action=hg_ghost --nparts=6"
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
ENDIF()


IF (${PROJECT_NAME}_ENABLE_SCOREC)
  IF (${PROJECT_NAME}_ENABLE_Pamgen)
    TRIBITS_ADD_TEST(
      pamgenMeshAdapterTest
      NAME pamgenMeshAdapterTest_ParMA
      NUM_MPI_PROCS 4
      COMM serial mpi
      ARGS
      "--action=parma"
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )
  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    APFMeshAdapterTest
    NAME APFMeshAdapterTest_ParMA
    SOURCES APFMeshAdapterTest.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "--action=parma --meshfile=4/ --modelfile=torus.dmg"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
  TRIBITS_ADD_TEST(
    APFMeshAdapterTest
    NAME APFMeshAdapterTest_ParMA_Vtx_Edge_Elm
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "--action=parma --parma_method=VtxEdgeElm --meshfile=4/ --modelfile=torus.dmg"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
  TRIBITS_ADD_TEST(
    APFMeshAdapterTest
    NAME APFMeshAdapterTest_ParMA_Ghost
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "--action=parma --parma_method=Ghost --meshfile=4/ --modelfile=torus.dmg"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
  TRIBITS_ADD_TEST(
    APFMeshAdapterTest
    NAME APFMeshAdapterTest_ParMA_Shape
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "--action=parma --parma_method=Shape --meshfile=4/ --modelfile=torus.dmg"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
  
  TRIBITS_ADD_TEST(
    APFMeshAdapterTest
    NAME APFMeshAdapterTest_ZOLTAN_HG
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "--action=zoltan_hg --meshfile=4/ --modelfile=torus.dmg"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
  TRIBITS_ADD_TEST(
    APFMeshAdapterTest
    NAME APFMeshAdapterTest_ZOLTAN_GHOST_HG
    NUM_MPI_PROCS 2
    COMM serial mpi
    ARGS
    "--action=hg_ghost --meshfile=pumiTri14/2/ --modelfile=pumiTri14/plate.dmg --layers=1"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    PartitionAndParMATest
    NAME ZoltanHGAndParMA_VtxElm
    SOURCES PartitionAndParMATest.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "--action=zoltan_hg --meshfile=4/ --modelfile=torus.dmg"
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

IF (Tpetra_INST_DOUBLE)
  # Double needed for file I/O.
  TRIBITS_ADD_EXECUTABLE(
    MapSimulater
    SOURCES TaskMappingSimulate.cpp
  )
ENDIF()

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
    torus.dmg
    4/0.smb
    4/1.smb	
    4/2.smb
    4/3.smb
    pumiTri14/plate.dmg
    pumiTri14/plate0.smb
    pumiTri14/2/0.smb
    pumiTri14/2/1.smb
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)









##########################################
IF (${PACKAGE_NAME}_ENABLE_Experimental)

IF (${PACKAGE_NAME}_ENABLE_Experimental_2D)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Partitioning2D
  SOURCES partition2DMatrix.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

ENDIF()

##########################################


##########################################

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PartitioningTree
  SOURCES partitioningTree.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

##########################################

IF (NOT Tpetra_ENABLE_CUDA)
IF (${PACKAGE_NAME}_ENABLE_Epetra)
IF (Tpetra_INST_INT_INT)

# Exercises the Epetra path for the XpetraMultiVector
# interfacing with the Multijagged algorithm.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    mj_epetra
    SOURCES mj_epetra.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

ENDIF()
ENDIF()
ENDIF()

