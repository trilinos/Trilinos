TRIBITS_ADD_EXECUTABLE(
solverfactory
SOURCES main.cpp
)

ADD_SUBDIRECTORIES(
ParameterLists
)

##############
# Amesos Klu #
##############
IF(HAVE_SHYLU_DDFROSCH_AMESOS)
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos_klu_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Amesos_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos_klu_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Amesos_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
ENDIF()
ENDIF()

###############
# Amesos2 Klu #
###############
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos2_klu_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Amesos2_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos2_klu_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Amesos2_Klu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos2_klu_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Amesos2_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos2_klu_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Amesos2_Klu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 1
)

###############
# Belos GMRES #
###############
IF(HAVE_SHYLU_DDFROSCH_BELOS)
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

#################
# Ifpack2 RILUK #
#################
IF(HAVE_SHYLU_DDFROSCH_IFPACK2)
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_ifpack2_riluk_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Ifpack2_RILUK.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_ifpack2_riluk_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Ifpack2_RILUK.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

####################################################
# FROSchPreconditioner TwoLevelBlockPreconditioner #
####################################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

###############################################
# FROSchPreconditioner TwoLevelPreconditioner #
###############################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

#########
# MueLu #
#########
IF(HAVE_SHYLU_DDFROSCH_MUELU)
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_muelu_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_MueLu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_muelu_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_MueLu.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_muelu_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_MueLu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_muelu_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_MueLu.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

#####################################################
# ThyraPreconditioner FROSch OneLevelPreconditioner #
#####################################################
IF(HAVE_SHYLU_DDFROSCH_THYRA AND HAVE_SHYLU_DDFROSCH_STRATIMIKOS)
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

###########################
# ThyraSolver Belos GMRES #
###########################
IF(HAVE_SHYLU_DDFROSCH_THYRA AND HAVE_SHYLU_DDFROSCH_STRATIMIKOS AND HAVE_SHYLU_DDFROSCH_BELOS)
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrasolver_belos_gmres_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraSolver_Belos_GMRES.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrasolver_belos_gmres_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraSolver_Belos_GMRES.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrasolver_belos_gmres_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraSolver_Belos_GMRES.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrasolver_belos_gmres_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraSolver_Belos_GMRES.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
