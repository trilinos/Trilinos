TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND
   ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MeshTying_cpOneDof
    SOURCE_FILES
      Q_mm.txt
      G_mm.txt
      GT_mm.txt
      C_mm.txt
      f_mm.txt
      Lagr2Dof.txt
      simple_1dof.xml
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MeshTying2D_cp_small
    SOURCE_FILES
      meshtying2D_small_matrix.mm
      meshtying2D_small_rhs.mm
      meshtying2D_small_nullspace1.mm
      meshtying2D_small_interface_dof_map_MPI1.mm
      meshtying2D_small_interface_dof_map_MPI4.mm
      meshtying2D_small_simple.xml
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MeshTying2D_cp_medium
    SOURCE_FILES
      meshtying2D_medium_matrix.mm
      meshtying2D_medium_rhs.mm
      meshtying2D_medium_nullspace1.mm
      meshtying2D_medium_interface_dof_map_MPI1.mm
      meshtying2D_medium_interface_dof_map_MPI4.mm
      meshtying2D_medium_simple.xml
      meshtying2D_medium_simple_spai.xml
    )

  TRIBITS_ADD_EXECUTABLE(
    MeshTyingBlocked_NodeBased
    SOURCES MeshTyingBlocked_NodeBased.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    MortarSurfaceCoupling_DofBased
    SOURCES MortarSurfaceCoupling_DofBased.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    MortarSurfaceCoupling_DofBased
    NAME "MortarSurfaceCoupling_DofBased_Blocked_SimpleSmoother_2dof_small"
    ARGS "--linAlgebra=Tpetra --xml=meshtying2D_small_simple.xml --nPrimalDofs=78 --nDualDofs=14 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_small --expectedNumIts=16"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    MortarSurfaceCoupling_DofBased
    NAME "MortarSurfaceCoupling_DofBased_Blocked_SimpleSmoother_2dof_medium"
    ARGS "--linAlgebra=Tpetra --xml=meshtying2D_medium_simple.xml --nPrimalDofs=530 --nDualDofs=54 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_medium --expectedNumIts=6"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    MortarSurfaceCoupling_DofBased
    NAME "MortarSurfaceCoupling_DofBased_Blocked_SimpleSmoother_2dof_medium"
    ARGS "--linAlgebra=Tpetra --xml=meshtying2D_medium_simple.xml --nPrimalDofs=530 --nDualDofs=54 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_medium --expectedNumIts=12"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    MortarSurfaceCoupling_DofBased
    NAME "MortarSurfaceCoupling_DofBased_Blocked_SimpleSmoother_2dof_medium_spai"
    ARGS "--linAlgebra=Tpetra --xml=meshtying2D_medium_simple_spai.xml --nPrimalDofs=530 --nDualDofs=54 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_medium --expectedNumIts=7"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    MeshTyingBlocked_NodeBased
    NAME "MeshTyingBlocked_NodeBased_simpleSmoother"
    ARGS "-linAlgebra=Tpetra"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

ENDIF()