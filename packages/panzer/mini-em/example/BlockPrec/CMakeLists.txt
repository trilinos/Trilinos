TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(BLOCK_PREC_SOURCES
  main.cpp
  MiniEM_helpers.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  BlockPrec
  SOURCES ${BLOCK_PREC_SOURCES}
  CATEGORIES BASIC PERFORMANCE
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBlockPrecFiles
  SOURCE_FILES
  maxwell.xml maxwell2D.xml
  maxwell-1stOrder.xml maxwell-pthOrder.xml
  maxwell-blob-R0.xml maxwell-blob-R1.xml maxwell-blob-R2.xml maxwell-blob-R3.xml maxwell-blob-R4.xml
  maxwell-bdot-small.xml maxwell-bdot-medium.xml maxwell-bdot-large.xml
  maxwell-exterior-small.xml maxwell-exterior-medium.xml maxwell-exterior-large.xml
  darcy.xml
  solverCG.xml solverGMRES.xml
  solverMueLuRefMaxwell.xml solverMueLuRefMaxwell2D.xml solverMueLuRefMaxwellEpetra.xml
  solverMueLuMaxwellHO.xml solverMueLuMaxwellHOCuda.xml
  solverMueLuRefMaxwellOpenMP.xml solverMueLuRefMaxwellCuda.xml
  solverAugmentationUseILU.xml
  solverMLRefMaxwell.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBlockPrecPerfFiles
  SOURCE_FILES
  maxwell-large.xml
  solverMueLuRefMaxwell.xml solverMueLuRefMaxwell2D.xml solverMueLuRefMaxwellEpetra.xml
  solverMueLuRefMaxwellOpenMP.xml solverMueLuRefMaxwellCuda.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  CATEGORIES PERFORMANCE
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

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_reuse"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=3 --linAlgebra=Tpetra --resetSolver"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_reuse"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=3 --linAlgebra=Tpetra --resetSolver"
   COMM mpi
   NUM_MPI_PROCS 4
   )

 # Large performance-testing runs of 3D BlockPrec_RefMaxwell (4 and 16 ranks)
TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Performance_4"
   ARGS "--stacked-timer --solver=MueLu-RefMaxwell --numTimeSteps=3 --linAlgebra=Tpetra --inputFile=maxwell-large.xml"
   COMM mpi
   NUM_MPI_PROCS 4
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Performance_16"
   ARGS "--stacked-timer --solver=MueLu-RefMaxwell --numTimeSteps=3 --linAlgebra=Tpetra --inputFile=maxwell-large.xml"
   COMM mpi
   NUM_MPI_PROCS 16
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Performance_16_reuse"
   ARGS "--stacked-timer --solver=MueLu-RefMaxwell --numTimeSteps=20 --linAlgebra=Tpetra --inputFile=maxwell-large.xml --resetSolver --test-name=\"MiniEM 3D RefMaxwell Reuse\""
   COMM mpi
   NUM_MPI_PROCS 16
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )


 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell2D"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=1 --linAlgebra=Tpetra --inputFile=maxwell2D.xml"
   COMM mpi
   NUM_MPI_PROCS 4
   )

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

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Epetra_reuse"
   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=3 --linAlgebra=Epetra --resetSolver"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )

 # TRIBITS_ADD_TEST(
 #   BlockPrec
 #   NAME "MiniEM-BlockPrec_RefMaxwell_Epetra_reuse"
 #   ARGS "--solver=MueLu-RefMaxwell --numTimeSteps=3 --linAlgebra=Epetra --resetSolver"
 #   COMM mpi
 #   NUM_MPI_PROCS 4
 #   )

ENDIF()

# CAG 5/16/22: Commented out tests that rely on 2nd order elements.
#              The specialization of 2nd order basis functions in
#              Intrepid2 is missing some functionality and will be
#              removed. These tests will work once 2nd order has been
#              replace by p-th order with p=2.
TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_MueLu_highOrder"
   ARGS
   # "--solver=MueLu-Maxwell-HO --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=2"
   "--solver=MueLu-Maxwell-HO --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3"
   # "--solver=MueLu-Maxwell-HO --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3 --pCoarsenSchedule=\"2,1\""
   # "--solver=MueLu-Maxwell-HO --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=2 --matrixFree"
   "--solver=MueLu-Maxwell-HO --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3 --matrixFree"
   # "--solver=MueLu-Maxwell-HO --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3 --pCoarsenSchedule=\"2,1\" --matrixFree"
   COMM serial mpi
   NUM_MPI_PROCS 4
   )

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


#################################################
# CG

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_CG_highOrder"
   ARGS
   "--solver=CG --linAlgebra=Tpetra --inputFile=maxwell-1stOrder.xml --basis-order=1"
   "--solver=CG --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=2"
   "--solver=CG --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3"
   COMM serial mpi
   NUM_MPI_PROCS 4
   )



#################################################
#################################################
# DARCY

#################################################
# CG

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_Darcy_CG"
   ARGS
   "--solver=CG --linAlgebra=Tpetra --inputFile=darcy.xml --basis-order=1 --x-elements=15  --y-elements=15 --z-elements=15"
   "--solver=CG --linAlgebra=Tpetra --inputFile=darcy.xml --basis-order=2 --x-elements=8  --y-elements=8 --z-elements=8"
   "--solver=CG --linAlgebra=Tpetra --inputFile=darcy.xml --basis-order=3 --x-elements=4  --y-elements=4 --z-elements=4"
   COMM serial mpi
   NUM_MPI_PROCS 4
   )
