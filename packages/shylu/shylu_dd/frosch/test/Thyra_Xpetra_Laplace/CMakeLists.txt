TRIBITS_ADD_EXECUTABLE(
thyra_xpetra_laplace
SOURCES main.cpp
)

ADD_SUBDIRECTORIES(
ParameterLists
)

######################
# GDSWPreconditioner #
######################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_GDSWP_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

#######################
# RGDSWPreconditioner #
#######################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_RGDSWP_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

###############################
# TwoLevelPreconditioner_GDSW #
###############################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_GDSW_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

################################
# TwoLevelPreconditioner_RGDSW #
################################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_RGDSW_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

#######################################
# TwoLevelPreconditioner_IPOUHarmonic #
#######################################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=1 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=1 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

####################################
# TwoLevelBlockPreconditioner_GDSW #
####################################

## Number of Blocks = 2

## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=2 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=2 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

## Number of Blocks = 3

## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD0_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD1_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN1_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN3_ORD0_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN1_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN3_ORD1_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD0_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=2 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD1_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --NB=3 --DPN=2 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN1_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=1 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN3_ORD0_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=3 --ORD=0 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN1_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=1 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_laplace
NAME test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM3_DPN3_ORD1_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --NB=3 --DPN=3 --ORD=1 --PLIST=ParameterLists/ParameterList_TwoLevelBlockPreconditioner_GDSW.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
