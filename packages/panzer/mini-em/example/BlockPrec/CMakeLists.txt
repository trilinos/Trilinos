INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(BLOCK_PREC_SOURCES
  main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  BlockPrec
  SOURCES ${BLOCK_PREC_SOURCES}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBlockPrecFiles
    SOURCE_FILES solverDefaultsAugmentation.xml solverDefaultsRefMaxwell.xml solverAugmentationUseILU.xml refMaxwellKokkos.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_Augmentation"
   ARGS "--numTimeSteps=1"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_Augmentation"
   ARGS "--numTimeSteps=1"
   COMM mpi
   NUM_MPI_PROCS 4
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell"
   ARGS "--numTimeSteps=1 --use-refmaxwell"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell"
   ARGS "--numTimeSteps=1 --use-refmaxwell"
   COMM mpi
   NUM_MPI_PROCS 4
   )
