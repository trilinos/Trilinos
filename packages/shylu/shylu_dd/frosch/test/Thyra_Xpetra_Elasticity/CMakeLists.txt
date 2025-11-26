TRIBITS_ADD_EXECUTABLE(
thyra_xpetra_elasticity
SOURCES main.cpp
CATEGORIES BASIC PERFORMANCE
)

ADD_SUBDIRECTORIES(
ParameterLists
)

######################
# GDSWPreconditioner #
######################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_GDSWP_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_GDSWP_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_GDSWPreconditioner_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

#######################
# RGDSWPreconditioner #
#######################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_RGDSWP_DIM2_TPETRA_SingleReduce
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu_Tpetra.xml"
COMM mpi
NUM_MPI_PROCS 4
)

IF (Trilinos_ENABLE_FLOAT)
    TRIBITS_ADD_TEST(
        thyra_xpetra_elasticity
        NAME test_thyra_xpetra_elasticity_RGDSWP_DIM2_TPETRA_Half
        ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_RGDSWPreconditioner_Klu_Half.xml"
        COMM mpi
        NUM_MPI_PROCS 4
    )
ENDIF ()

###############################
# TwoLevelPreconditioner_GDSW #
###############################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_GDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_GDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_GDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

################################
# TwoLevelPreconditioner_RGDSW #
################################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_RGDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_RGDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

############################################
# TwoLevelPreconditioner_IPOUHarmonic_GDSW #
############################################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

IF(HAVE_SHYLU_DDFROSCH_ZOLTAN2)
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_ML_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu_ML.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_ML_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSW_Klu_ML.xml "
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

################################################
# TwoLevelPreconditioner_IPOUHarmonic_GDSWStar #
################################################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_GDSWStar_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

#############################################
# TwoLevelPreconditioner_IPOUHarmonic_RGDSW #
#############################################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM3_TPETRA
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 8
)

IF(HAVE_SHYLU_DDFROSCH_ZOLTAN2)
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_ML_DIM2_TPETRA
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_ML.xml"
COMM mpi
NUM_MPI_PROCS 4
)
ENDIF()


### Performance Tests ###
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME performance_TLP_IPOUHarmonic_RGDSW_DIM2
ARGS "--M=100 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 4
CATEGORIES PERFORMANCE
RUN_SERIAL
)
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME performance_TLP_IPOUHarmonic_RGDSW_DIM2
ARGS "--M=100 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 9
CATEGORIES PERFORMANCE
RUN_SERIAL
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME performance_TLP_IPOUHarmonic_RGDSW_DIM2
ARGS "--M=100 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 16
CATEGORIES PERFORMANCE
RUN_SERIAL
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME performance_TLP_IPOUHarmonic_RGDSW_DIM2
ARGS "--M=100 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu.xml"
COMM mpi
NUM_MPI_PROCS 25
CATEGORIES PERFORMANCE
RUN_SERIAL
)



## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_TPETRA_DropCoupling
ARGS "--M=3 --DIM=2 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_2D.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_elasticity
NAME test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM3_TPETRA_DropCoupling
ARGS "--M=3 --DIM=3 --O=1 --PLIST=ParameterLists/ParameterList_TwoLevelPreconditioner_IPOUHarmonic_RGDSW_Klu_Drop_Coupling_3D.xml"
COMM mpi
NUM_MPI_PROCS 8
)
