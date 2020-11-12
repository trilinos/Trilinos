TRIBITS_ADD_EXECUTABLE(
thyra_xpetra_elasticity
SOURCES main.cpp
)

ADD_SUBDIRECTORIES(
ParameterLists
)

######################
# GDSWPreconditioner #
######################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_GDSWP_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_GDSWP_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_GDSWP_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_GDSWP_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

#######################
# RGDSWPreconditioner #
#######################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

###############################
# TwoLevelPreconditioner_GDSW #
###############################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_GDSW_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_GDSW_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_GDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_GDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

################################
# TwoLevelPreconditioner_RGDSW #
################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_RGDSW_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_RGDSW_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_RGDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_RGDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

############################################
# TwoLevelPreconditioner_IPOUHarmonic_GDSW #
############################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

IF(HAVE_SHYLU_DDFROSCH_ZOLTAN2)
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_ML_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu_ML.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_ML_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu_ML.xml  --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

################################################
# TwoLevelPreconditioner_IPOUHarmonic_GDSWStar #
################################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

#############################################
# TwoLevelPreconditioner_IPOUHarmonic_RGDSW #
#############################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_EPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM3_EPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

IF(HAVE_SHYLU_DDFROSCH_ZOLTAN2)
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_ML_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_ML.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
ENDIF()

## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_EPETRA_DropCoupling
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_2D.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM3_EPETRA_DropCoupling
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_3D.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_TPETRA_DropCoupling
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_2D.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM3_TPETRA_DropCoupling
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_3D.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
