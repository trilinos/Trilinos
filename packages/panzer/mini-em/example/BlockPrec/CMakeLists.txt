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
  maxwell.xml maxwell2D.xml
  maxwell-blob-R0.xml maxwell-blob-R1.xml maxwell-blob-R2.xml maxwell-blob-R3.xml maxwell-blob-R4.xml
  maxwell-bdot-small.xml maxwell-bdot-medium.xml maxwell-bdot-large.xml
  maxwell-exterior-small.xml maxwell-exterior-medium.xml maxwell-exterior-large.xml
  solverAugmentation.xml solverAugmentationEpetra.xml
  solverCG.xml
  solverMueLuRefMaxwell.xml solverMueLuRefMaxwell2D.xml solverMueLuRefMaxwellEpetra.xml
  solverMueLuRefMaxwellOpenMP.xml solverMueLuRefMaxwellCuda.xml
  solverAugmentationUseILU.xml
  solverMLRefMaxwell.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

# rather than list all files like above, I glob all xmls
FILE(GLOB ${PACKAGE_NAME}_InputFiles *.xml)

# allow the user to override the install location and
# provide a simple defintion for the default
# EXAMPLE_INSTALL_DIR/Panzer/
IF (NOT DEFINED ${PACKAGE_NAME}_INSTALL_PREFIX)
  SET(${PACKAGE_NAME}_INSTALL_PREFIX "${${PROJECT_NAME}_INSTALL_EXAMPLE_DIR}/${PACKAGE_NAME}/")
ELSE()
  # inform the user where it will go, since they overrode the default
  PRINT_VAR(${PACKAGE_NAME}_INSTALL_PREFIX)
ENDIF()

INSTALL(FILES ${${PACKAGE_NAME}_InputFiles}
        DESTINATION ${${PACKAGE_NAME}_INSTALL_PREFIX}
       )
INSTALL(TARGETS "${PACKAGE_NAME}_BlockPrec"
        COMPONENT ${PACKAGE_NAME} RUNTIME DESTINATION ${${PACKAGE_NAME}_INSTALL_PREFIX})

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

IF(NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_OPENMP AND MueLu_ENABLE_Epetra AND NOT HAVE_TPETRA_INST_INT_LONG_LONG AND HAVE_TPETRA_INST_INT_INT)
 # We want GlobalOrdinal=int. Panzer picks long long, if it is available, and otherwise Tpetra's default.

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

IF(NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_OPENMP)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_RefMaxwell2D"
    ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Tpetra --inputFile=maxwell2D.xml"
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()

IF(NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_OPENMP AND MueLu_ENABLE_Epetra AND NOT HAVE_TPETRA_INST_INT_LONG_LONG AND HAVE_TPETRA_INST_INT_INT)
 # We want GlobalOrdinal=int. Panzer picks long long, if it is available, and otherwise Tpetra's default.

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Epetra"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Epetra"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Epetra"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Epetra"
   COMM mpi
   NUM_MPI_PROCS 4
   )

ENDIF()

#################################################
# ML RefMaxwell solver

IF (NOT Kokkos_ENABLE_CUDA AND ${PACKAGE_NAME}_ENABLE_ML AND MueLu_ENABLE_Epetra)

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
