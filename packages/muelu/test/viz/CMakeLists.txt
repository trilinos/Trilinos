# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
)

#
# Executable
#

#TODO: how to remove that?
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

#
# Tests
#

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    AggExport
    SOURCES AggExport.cpp
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    AggExport
    NAME "VizElasticity3DPointCloudTpetra"
    ARGS "--matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --linAlgebra=Tpetra --xml=vizElasticity3DPointCloud.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    AggExport
    NAME "VizLaplace3DPointCloudTpetra"
    ARGS "--matrixType=Laplace3D --nx=10 --ny=10 --nz=10 --linAlgebra=Tpetra --xml=vizLaplace3DPointCloud.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    AggExport
    NAME "VizLaplace3DPointCloudAggregateQualityTpetra"
    ARGS "--matrixType=Laplace3D --nx=10 --ny=10 --nz=10 --linAlgebra=Tpetra --xml=vizLaplace3DPointCloudAggregateQuality.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    AggExport
    NAME "VizElasticity3DConvexHullsTpetra"
    ARGS "--matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --linAlgebra=Tpetra --xml=vizElasticity3DConvexHulls.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    AggExport
    NAME "VizLaplace3DConvexHullsTpetra"
    ARGS "--matrixType=Laplace3D --nx=10 --ny=10 --nz=10 --linAlgebra=Tpetra --xml=vizLaplace3DConvexHulls.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    AggExport
    NAME "VizLaplace3DConvexHullsAggregateQualityTpetra"
    ARGS "--matrixType=Laplace3D --nx=10 --ny=10 --nz=10 --linAlgebra=Tpetra --xml=vizLaplace3DConvexHullsAggregateQuality.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  IF(${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_AztecOO)

    TRIBITS_ADD_EXECUTABLE(
      Viz
      SOURCES Viz.cpp
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Viz
      NAME "Viz3DTpetra"
      ARGS "--matrixType=Laplace3D --nx=20 --ny=20 --nz=20 --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Viz
      NAME "Viz2DTpetra"
      ARGS "--matrixType=Laplace2D --nx=40 --ny=40 --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
    )

    # Test Epetra variant, too.
    IF(${PACKAGE_NAME}_ENABLE_Epetra
       AND ${PACKAGE_NAME}_ENABLE_EpetraExt
       AND ${PACKAGE_NAME}_ENABLE_Ifpack
       AND ${PACKAGE_NAME}_ENABLE_Amesos
       AND (NOT Xpetra_INT_LONG_LONG)
    )

      MUELU_ADD_SERIAL_AND_MPI_TEST(
        Viz
        NAME "Viz3DEpetra"
        ARGS "--matrixType=Laplace3D --nx=20 --ny=20 --nz=20 --linAlgebra=Epetra"
        NUM_MPI_PROCS 4
        COMM serial mpi
      )

      MUELU_ADD_SERIAL_AND_MPI_TEST(
        Viz
        NAME "Viz2DEpetra"
        ARGS "--matrixType=Laplace2D --nx=40 --ny=40 --linAlgebra=Epetra"
        NUM_MPI_PROCS 4
        COMM serial mpi
      )

    ENDIF() #endif Epetra

  ENDIF()

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Viz_cp
  SOURCE_FILES
    Output/MPI-Viz-Output-Elasticity3D-ConvexHulls-Level0-Proc0.gold
    Output/MPI-Viz-Output-Elasticity3D-ConvexHulls-Level0-Proc1.gold
    Output/MPI-Viz-Output-Elasticity3D-ConvexHulls-Level0-Proc2.gold
    Output/MPI-Viz-Output-Elasticity3D-ConvexHulls-Level0-Proc3.gold
    Output/MPI-Viz-Output-Elasticity3D-ConvexHulls-Level0.gold
    Output/MPI-Viz-Output-Elasticity3D-PointCloud-Level0-Proc0.gold
    Output/MPI-Viz-Output-Elasticity3D-PointCloud-Level0-Proc1.gold
    Output/MPI-Viz-Output-Elasticity3D-PointCloud-Level0-Proc2.gold
    Output/MPI-Viz-Output-Elasticity3D-PointCloud-Level0-Proc3.gold
    Output/MPI-Viz-Output-Elasticity3D-PointCloud-Level0.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-AggregateQuality-Level0-Proc0.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-AggregateQuality-Level0-Proc1.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-AggregateQuality-Level0-Proc2.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-AggregateQuality-Level0-Proc3.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-AggregateQuality-Level0.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-Level0-Proc0.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-Level0-Proc1.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-Level0-Proc2.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-Level0-Proc3.gold
    Output/MPI-Viz-Output-Laplace3D-ConvexHulls-Level0.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-AggregateQuality-Level0-Proc0.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-AggregateQuality-Level0-Proc1.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-AggregateQuality-Level0-Proc2.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-AggregateQuality-Level0-Proc3.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-AggregateQuality-Level0.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-Level0-Proc0.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-Level0-Proc1.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-Level0-Proc2.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-Level0-Proc3.gold
    Output/MPI-Viz-Output-Laplace3D-PointCloud-Level0.gold
    vizElasticity3DConvexHulls.xml
    vizElasticity3DPointCloud.xml
    vizLaplace3DConvexHulls.xml
    vizLaplace3DConvexHullsAggregateQuality.xml
    vizLaplace3DPointCloud.xml
    vizLaplace3DPointCloudAggregateQuality.xml
    viztest.xml
)
