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
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos2_klu_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Amesos2_Klu.xml"
COMM mpi
NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_amesos2_klu_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Amesos2_Klu.xml"
COMM mpi
NUM_MPI_PROCS 1
)

###############
# Belos GMRES #
###############
IF(HAVE_SHYLU_DDFROSCH_BELOS)
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_belos_gmres_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Belos_GMRES.xml"
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
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_Ifpack2_RILUK.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_ifpack2_riluk_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_Ifpack2_RILUK.xml"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

####################################################
# FROSchPreconditioner TwoLevelBlockPreconditioner #
####################################################
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelBlockPreconditioner.xml"
COMM mpi
NUM_MPI_PROCS 8
)

###############################################
# FROSchPreconditioner TwoLevelPreconditioner #
###############################################
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_FROSchPreconditioner_TwoLevelPreconditioner.xml"
COMM mpi
NUM_MPI_PROCS 8
)

#########
# MueLu #
#########
IF(HAVE_SHYLU_DDFROSCH_MUELU)
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_muelu_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_MueLu.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_muelu_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_MueLu.xml"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

#####################################################
# ThyraPreconditioner FROSch OneLevelPreconditioner #
#####################################################
IF(HAVE_SHYLU_DDFROSCH_THYRA AND HAVE_SHYLU_DDFROSCH_STRATIMIKOS)
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraPreconditioner_FROSch_OneLevelPreconditioner.xml"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

###########################
# ThyraSolver Belos GMRES #
###########################
IF(HAVE_SHYLU_DDFROSCH_THYRA AND HAVE_SHYLU_DDFROSCH_STRATIMIKOS AND HAVE_SHYLU_DDFROSCH_BELOS)
##Tpetra
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrasolver_belos_gmres_DIM2_TPETRA
ARGS "--DIM=2 --PLIST=ParameterLists/ParameterList_ThyraSolver_Belos_GMRES.xml"
COMM mpi
NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
solverfactory
NAME test_solverfactory_thyrasolver_belos_gmres_DIM3_TPETRA
ARGS "--DIM=3 --PLIST=ParameterLists/ParameterList_ThyraSolver_Belos_GMRES.xml"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()
