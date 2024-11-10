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
  maxwell-donut.xml maxwell-analyticSolution.xml
  darcyTet.xml darcyHex.xml
  darcyTetAnalytic.xml darcyHexAnalytic.xml
  solverCG.xml solverGMRES.xml
  solverMueLu.xml solverMueLu2D.xml solverMueLuEpetra.xml
  solverMueLuOpenMP.xml solverMueLuCuda.xml solverMueLuTruncated.xml solverMueLuTPL.xml
  solverMueLuHO.xml solverMueLuHOCuda.xml
  solverAugmentationUseILU.xml
  solverML.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBlockPrecPerfFiles
  SOURCE_FILES
  maxwell-large.xml
  solverMueLu.xml solverMueLu2D.xml solverMueLuEpetra.xml
  solverMueLuOpenMP.xml solverMueLuCuda.xml
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
#################################################
# MAXELL

#################################################
# MueLu RefMaxwell solver

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Maxwell_MueLu"
   POSTFIX_AND_ARGS_0 "order1"       --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra
   POSTFIX_AND_ARGS_1 "order1_reuse" --solver=MueLu --numTimeSteps=3 --linAlgebra=Tpetra --resetSolver
   NUM_MPI_PROCS 1
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Maxwell_MueLu"
   POSTFIX_AND_ARGS_0 "order1"              --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra
   POSTFIX_AND_ARGS_1 "order1_reuse"        --solver=MueLu --numTimeSteps=3 --linAlgebra=Tpetra --resetSolver
   POSTFIX_AND_ARGS_2 "order2"              --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=2 --pCoarsenSchedule="1"
   POSTFIX_AND_ARGS_3 "order3"              --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3 --pCoarsenSchedule="1"
   POSTFIX_AND_ARGS_4 "order3,2"            --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3 --pCoarsenSchedule="2,1"
   POSTFIX_AND_ARGS_5 "order2_matrixFree"   --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=2  --pCoarsenSchedule="1" --matrixFree
   POSTFIX_AND_ARGS_6 "order3_matrixFree"   --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3  --pCoarsenSchedule="1" --matrixFree
   POSTFIX_AND_ARGS_7 "order3,2_matrixFree" --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3 --pCoarsenSchedule="2,1" --matrixFree
   POSTFIX_AND_ARGS_8 "2D"                  --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra --inputFile=maxwell2D.xml
   POSTFIX_AND_ARGS_9 "order1_analytic"     --solver=MueLu --linAlgebra=Tpetra --inputFile=maxwell-analyticSolution.xml
   POSTFIX_AND_ARGS_10 "order1_truncated"    --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra --truncateMueLuHierarchy
   NUM_MPI_PROCS 4
   )


IF (${PACKAGE_NAME}_ENABLE_ShyLU_NodeTacho)
  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "Maxwell_MueLu"
    POSTFIX_AND_ARGS_0 "order1_tpl"    --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra --tpl
    ENVIRONMENT TPETRA_FUSED_RESIDUAL=OFF
    NUM_MPI_PROCS 4
  )
 ENDIF()


 # Large performance-testing runs of 3D BlockPrec_RefMaxwell (4 and 16 ranks)
TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Performance_4"
   ARGS "--stacked-timer --solver=MueLu --numTimeSteps=3 --linAlgebra=Tpetra --inputFile=maxwell-large.xml"
   COMM mpi
   NUM_MPI_PROCS 4
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Performance_16"
   ARGS "--stacked-timer --solver=MueLu --numTimeSteps=3 --linAlgebra=Tpetra --inputFile=maxwell-large.xml"
   COMM mpi
   NUM_MPI_PROCS 16
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )

IF (${PACKAGE_NAME}_ENABLE_ShyLU_NodeTacho)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_RefMaxwell_Performance_4_tpl"
    ARGS "--stacked-timer --solver=MueLu --numTimeSteps=3 --linAlgebra=Tpetra --inputFile=maxwell-large.xml --tpl --test-name=\"MiniEM 3D RefMaxwell TPL\""
    ENVIRONMENT TPETRA_FUSED_RESIDUAL=OFF
    COMM mpi
    NUM_MPI_PROCS 4
    CATEGORIES PERFORMANCE
    RUN_SERIAL
  )

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "MiniEM-BlockPrec_RefMaxwell_Performance_16_tpl"
    ARGS "--stacked-timer --solver=MueLu --numTimeSteps=3 --linAlgebra=Tpetra --inputFile=maxwell-large.xml --tpl --test-name=\"MiniEM 3D RefMaxwell TPL\""
    ENVIRONMENT TPETRA_FUSED_RESIDUAL=OFF
    COMM mpi
    NUM_MPI_PROCS 16
    CATEGORIES PERFORMANCE
    RUN_SERIAL
  )

ENDIF()

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "MiniEM-BlockPrec_RefMaxwell_Performance_16_reuse"
   ARGS "--stacked-timer --solver=MueLu --numTimeSteps=20 --linAlgebra=Tpetra --inputFile=maxwell-large.xml --resetSolver --test-name=\"MiniEM 3D RefMaxwell Reuse\""
   COMM mpi
   NUM_MPI_PROCS 16
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )


IF(NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_OPENMP AND MueLu_ENABLE_Epetra AND NOT HAVE_TPETRA_INST_INT_LONG_LONG AND HAVE_TPETRA_INST_INT_INT)
 # We want GlobalOrdinal=int. Panzer picks long long, if it is available, and otherwise Tpetra's default.

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Maxwell_MueLu_Epetra"
   POSTFIX_AND_ARGS_0 "order1"       --solver=MueLu --numTimeSteps=1 --linAlgebra=Epetra
   POSTFIX_AND_ARGS_1 "order1_reuse" --solver=MueLu --numTimeSteps=3 --linAlgebra=Epetra --resetSolver
   NUM_MPI_PROCS 1
   )

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Maxwell_MueLu_Epetra"
   ARGS "--solver=MueLu --numTimeSteps=1 --linAlgebra=Epetra"
   NUM_MPI_PROCS 4
   )

 # TRIBITS_ADD_TEST(
 #   BlockPrec
 #   NAME "RefMaxwell_Epetra_reuse"
 #   ARGS "--solver=MueLu --numTimeSteps=3 --linAlgebra=Epetra --resetSolver"
 #   COMM mpi
 #   NUM_MPI_PROCS 4
 #   )

ENDIF()


#################################################
# ML RefMaxwell solver

IF (NOT Kokkos_ENABLE_CUDA AND ${PACKAGE_NAME}_ENABLE_ML AND MueLu_ENABLE_Epetra)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "Maxwell_ML"
    ARGS "--solver=ML --numTimeSteps=1 --linAlgebra=Epetra"
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "Maxwell_ML"
    ARGS "--solver=ML --numTimeSteps=1 --linAlgebra=Epetra"
    NUM_MPI_PROCS 4
    )

ENDIF()


#################################################
# CG

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Maxwell_CG"
   POSTFIX_AND_ARGS_0 "order1" --solver=CG --linAlgebra=Tpetra --inputFile=maxwell-1stOrder.xml --basis-order=1
   POSTFIX_AND_ARGS_1 "order2" --solver=CG --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=2
   POSTFIX_AND_ARGS_2 "order3" --solver=CG --linAlgebra=Tpetra --inputFile=maxwell-pthOrder.xml --basis-order=3
   NUM_MPI_PROCS 4
   )



#################################################
#################################################
# DARCY

#################################################
# MueLu graddiv solver

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Darcy_MueLu"
   POSTFIX_AND_ARGS_0
   "order1" --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra --inputFile=darcyTet.xml --x-elements=15  --y-elements=15 --z-elements=15
   NUM_MPI_PROCS 1
   RUN_SERIAL
   )

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Darcy_MueLu"
   POSTFIX_AND_ARGS_0 "order1_tet"     --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra --inputFile=darcyTet.xml --x-elements=15  --y-elements=15 --z-elements=15
   POSTFIX_AND_ARGS_1 "order1_hex"     --solver=MueLu --numTimeSteps=1 --linAlgebra=Tpetra --inputFile=darcyHex.xml --x-elements=15  --y-elements=15 --z-elements=15
   POSTFIX_AND_ARGS_2 "order2"     --solver=MueLu --linAlgebra=Tpetra --inputFile=darcyTet.xml --basis-order=2 --pCoarsenSchedule=\"1\" --x-elements=8  --y-elements=8 --z-elements=8
   POSTFIX_AND_ARGS_3 "order3"         --solver=MueLu --linAlgebra=Tpetra --inputFile=darcyTet.xml --basis-order=3 --pCoarsenSchedule=\"1\" --x-elements=4  --y-elements=4 --z-elements=4
   POSTFIX_AND_ARGS_4 "order3,2"   --solver=MueLu --linAlgebra=Tpetra --inputFile=darcyTet.xml --basis-order=3 --pCoarsenSchedule=\"2,1\" --x-elements=4  --y-elements=4 --z-elements=4
   NUM_MPI_PROCS 4
   )


IF(NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_OPENMP AND MueLu_ENABLE_Epetra AND NOT HAVE_TPETRA_INST_INT_LONG_LONG AND HAVE_TPETRA_INST_INT_INT)
 # We want GlobalOrdinal=int. Panzer picks long long, if it is available, and otherwise Tpetra's default.

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Darcy_MueLu_Epetra_order1"
   ARGS "--solver=MueLu --numTimeSteps=1 --linAlgebra=Epetra --inputFile=darcyTet.xml"
   NUM_MPI_PROCS 1
 )

 TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Darcy_MueLu_Epetra_order1"
   ARGS "--solver=MueLu --numTimeSteps=1 --linAlgebra=Epetra --inputFile=darcyTet.xml"
   NUM_MPI_PROCS 4
 )

ENDIF()

#################################################
# ML graddiv solver

IF (NOT Kokkos_ENABLE_CUDA AND ${PACKAGE_NAME}_ENABLE_ML AND MueLu_ENABLE_Epetra)

  TRIBITS_ADD_TEST(
    BlockPrec
    NAME "Darcy_ML"
    ARGS "--solver=ML --numTimeSteps=1 --linAlgebra=Epetra --inputFile=darcyHex.xml  --x-elements=15  --y-elements=15 --z-elements=15"
    NUM_MPI_PROCS 1
    )

  # TRIBITS_ADD_TEST(
  #   BlockPrec
  #   NAME "Darcy_ML"
  #   ARGS "--solver=ML --numTimeSteps=1 --linAlgebra=Epetra --inputFile=darcyHex.xml  --x-elements=15  --y-elements=15 --z-elements=15"
  #   NUM_MPI_PROCS 4
  #   )

ENDIF()


#################################################
# CG

TRIBITS_ADD_TEST(
   BlockPrec
   NAME "Darcy_CG"
   POSTFIX_AND_ARGS_0 "order1" --solver=CG --linAlgebra=Tpetra --inputFile=darcyTet.xml --basis-order=1 --x-elements=15 --y-elements=15 --z-elements=15
   POSTFIX_AND_ARGS_1 "order2" --solver=CG --linAlgebra=Tpetra --inputFile=darcyTet.xml --basis-order=2 --x-elements=8  --y-elements=8  --z-elements=8
   POSTFIX_AND_ARGS_2 "order3" --solver=CG --linAlgebra=Tpetra --inputFile=darcyTet.xml --basis-order=3 --x-elements=4  --y-elements=4  --z-elements=4
   NUM_MPI_PROCS 4
   )
