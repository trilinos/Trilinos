TRIBITS_ADD_EXECUTABLE(
solverfactory
SOURCES main.cpp
)

ADD_SUBDIRECTORIES(
ParameterLists
)

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
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 1
)

#####################################################
# ThyraPreconditioner FROSch OneLevelPreconditioner #
#####################################################
## Epetra
IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditionerDIM2_EPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditionerDIM3_EPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditionerDIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditionerDIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)


###########################
# ThyraSolver Belos GMRES #
###########################
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
