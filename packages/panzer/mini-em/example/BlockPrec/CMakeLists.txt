INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(BLOCK_PREC_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  BlockPrec
  SOURCES ${BLOCK_PREC_SOURCES}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBlockPrecFiles
  SOURCE_FILES
  solverAugmentation.xml solverAugmentationEpetra.xml
  solverMueLuRefMaxwell.xml solverMueLuRefMaxwellEpetra.xml
  solverMueLuRefMaxwellOpenMP.xml solverMueLuRefMaxwellCuda.xml
  solverAugmentationUseILU.xml
  solverMLRefMaxwell.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

#################################################
# Augmentation solver

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_Augmentation"
   ARGS "--solver=Augmentation --numTimeSteps=1 --linAlgebra=Tpetra"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_Augmentation"
   ARGS "--solver=Augmentation --numTimeSteps=1 --linAlgebra=Tpetra"
   COMM mpi
   NUM_MPI_PROCS 4
   )

IF(NOT Kokkos_ENABLE_Cuda AND NOT Kokkos_ENABLE_OpenMP)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_Augmentation_Epetra"
    ARGS "--solver=Augmentation --numTimeSteps=1 --linAlgebra=Epetra"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_Augmentation_Epetra"
    ARGS "--solver=Augmentation --numTimeSteps=1 --linAlgebra=Epetra"
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()


#################################################
# MueLu RefMaxwell solver

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Tpetra"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Tpetra"
   COMM mpi
   NUM_MPI_PROCS 4
   )

IF(NOT Kokkos_ENABLE_Cuda AND NOT Kokkos_ENABLE_OpenMP)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_RefMaxwell_Epetra"
    ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Epetra"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  # TRIBITS_ADD_TEST(
  #   BlockPrec
  #   NAME "MiniEM-BlockPrec_RefMaxwell_Epetra"
  #   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Epetra"
  #   COMM mpi
  #   NUM_MPI_PROCS 4
  #   )

ENDIF()

#################################################
# ML RefMaxwell solver

IF (NOT Kokkos_ENABLE_Cuda AND ${PACKAGE_NAME}_ENABLE_ML)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_ML-RefMaxwell"
    ARGS "--solver=ML-RefMaxwell --numTimeSteps=1 --linAlgebra=Epetra"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_ML-RefMaxwell"
    ARGS "--solver=ML-RefMaxwell --numTimeSteps=1 --linAlgebra=Epetra"
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()
