
#
# Tests for src/algorithms
#   For now - we test algorithms with system tests.

#
# Tests for src/environment
#
#
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Machine
  SOURCES environment/Machine.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TimerManager
  SOURCES environment/TimerManager.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DebugManager
  SOURCES environment/DebugManager.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Environment
  SOURCES environment/Environment.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MetricOutputManager
  SOURCES environment/MetricOutputManager.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Parameters
  SOURCES environment/Parameters.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AllParameters
  SOURCES environment/AllParameters.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InputTraitsGood
  SOURCES input/InputTraitsGood.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InputTraitsBad
  SOURCES input/InputTraitsBad.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
#
# Tests for src/input
#
#
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicCoordinateInput
  SOURCES input/BasicCoordinateInput.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
 
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicIdentifierInput
  SOURCES input/BasicIdentifierInput.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicKokkosIdentifierInput
  SOURCES input/BasicKokkosIdentifierInput.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicVectorInput
  SOURCES input/BasicVectorInput.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PROJECT_NAME}_ENABLE_Pamgen)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    PamgenMeshInput
    SOURCES input/PamgenMeshInput.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

IF (${PROJECT_NAME}_ENABLE_SCOREC)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    APFMeshInput
    SOURCES input/APFMeshInput.cpp
    NUM_MPI_PROCS 2
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_scorec_unit_tests
  SOURCE_FILES
    pumiTri14/2/0.smb
    pumiTri14/2/1.smb
    pumiTri14/plate.dmg
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../partition
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

ENDIF()
 
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraRowMatrixInput
  SOURCES input/TpetraRowMatrixInput.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraRowGraphInput
  SOURCES input/TpetraRowGraphInput.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraTraits
  SOURCES input/XpetraTraits.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
 
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraCrsMatrixInput
  SOURCES input/XpetraCrsMatrixInput.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
 
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraCrsGraphInput
  SOURCES input/XpetraCrsGraphInput.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraVectorInput
  SOURCES input/XpetraVectorInput.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE(
  XpetraMultiVectorInput
  SOURCES input/XpetraMultiVectorInput.cpp 
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  XpetraMultiVectorInput
  NAME XpetraMultiVectorInput_noweights
  ARGS "--file=simple"
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  XpetraMultiVectorInput
  NAME XpetraMultiVectorInput_oneweight
  ARGS "--file=vwgt"
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  XpetraMultiVectorInput
  NAME XpetraMultiVectorInput_twoweights
  ARGS "--file=vwgt2"
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

#
# Tests for src/models
#
#
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GraphModel
  SOURCES models/GraphModel.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PROJECT_NAME}_ENABLE_SCOREC)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    HyperGraphModel
    SOURCES models/HyperGraphModel.cpp
    NUM_MPI_PROCS 2
    COMM mpi serial
  )
ENDIF()

#
# Tests for src/models
#
#
IF (${PROJECT_NAME}_ENABLE_Pamgen)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    GraphModel2ndAdjsFromAdjs
    SOURCES models/GraphModel2ndAdjsFromAdjs.cpp 
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IdentifierModel
  SOURCES models/IdentifierModel.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CoordinateModel
  SOURCES models/CoordinateModel.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

#
# Tests for src/problems
#
#    We test the Problem objects with system tests, and the
#    Solution objects with unit tests.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PartitioningSolution
  SOURCES problems/PartitioningSolution.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Mapping
  SOURCES problems/Mapping.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

#
# Tests for src/util
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StridedData
  SOURCES util/StridedData.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Metric
  SOURCES util/Metric.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TPLTraits
  SOURCES util/TPLTraits.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AlltoAll
  SOURCES util/AlltoAll.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  componentMetrics
  SOURCES util/componentMetrics.cpp
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  findUniqueGids
  SOURCES util/findUniqueGids.cpp
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_unit_tests
  SOURCE_FILES
    Poisson.xml
    simple.mtx
    simple_coord.mtx
    simple.graph
    simple.coords
    vwgt.graph
    vwgt.coords
    vwgt2.graph
    vwgt2.coords
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
