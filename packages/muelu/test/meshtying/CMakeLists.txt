INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND
    ${PACKAGE_NAME}_ENABLE_Ifpack2 AND
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
      meshtying2D_small_interface_dof_map.mm
      meshtying2D_small_simple.xml
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MeshTying2D_cp_medium
    SOURCE_FILES
      meshtying2D_medium_matrix.mm
      meshtying2D_medium_rhs.mm
      meshtying2D_medium_nullspace1.mm
      meshtying2D_medium_interface_dof_map.mm
      meshtying2D_medium_simple.xml
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MeshTyingBlocked_SimpleSmoother
    SOURCES MeshTyingBlocked_SimpleSmoother.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    MortarSurfaceCoupling_DofBased
    SOURCES MortarSurfaceCoupling_DofBased.cpp
    COMM serial mpi
    )

  # The MortarSurfaceCoupling tests only work with real Scalar types,
  # so choose an instantiation with real scalars (if any is available)
  SET(MortarSurfaceCoupling_Instatiation "")
  IF(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
    SET(MortarSurfaceCoupling_Instatiation "DOUBLE_INT_LONGLONGINT")
  ELSEIF(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
    SET(MortarSurfaceCoupling_Instatiation "DOUBLE_INT_LONGINT")
  ELSEIF(HAVE_MUELU_INST_DOUBLE_INT_INT)
    SET(MortarSurfaceCoupling_Instatiation "DOUBLE_INT_INT")
  ELSEIF(HAVE_MUELU_INST_FLOAT_INT_INT)
    SET(MortarSurfaceCoupling_Instatiation "FLOAT_INT_INT")
  ENDIF()

  IF(MortarSurfaceCoupling_Instatiation)
    TRIBITS_ADD_TEST(
      MortarSurfaceCoupling_DofBased
      NAME "MeshTyingBlocked_SimpleSmoother_2dof_small"
      ARGS "--linAlgebra=Tpetra --xml=meshtying2D_small_simple.xml --nPrimalDofs=78 --nDualDofs=14 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_small --expectedNumIts=16 --instantiation=${MortarSurfaceCoupling_Instatiation}"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_TEST(
      MortarSurfaceCoupling_DofBased
      NAME "MeshTyingBlocked_SimpleSmoother_2dof_medium"
      ARGS "--linAlgebra=Tpetra --xml=meshtying2D_medium_simple.xml --nPrimalDofs=530 --nDualDofs=54 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_medium --expectedNumIts=6 --instantiation=${MortarSurfaceCoupling_Instatiation}"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_TEST(
      MortarSurfaceCoupling_DofBased
      NAME "MeshTyingBlocked_SimpleSmoother_2dof_medium"
      ARGS "--linAlgebra=Tpetra --xml=meshtying2D_medium_simple.xml --nPrimalDofs=530 --nDualDofs=54 --numPrimalDofsPerNode=2 --numDualDofsPerNode=2 --probName=meshtying2D_medium --expectedNumIts=11 --instantiation=${MortarSurfaceCoupling_Instatiation}"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
  ENDIF()

ENDIF()
