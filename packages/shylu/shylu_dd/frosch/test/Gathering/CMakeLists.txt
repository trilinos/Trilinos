TRIBITS_ADD_EXECUTABLE(
gathering
SOURCES main.cpp
)

ADD_SUBDIRECTORIES(
ParameterLists
)

######################
# GDSWCoarseOperator #
######################
## Tpetra
TRIBITS_ADD_TEST(
gathering
NAME test_gathering_TPETRA
ARGS "--M=4 --DIM=3 --ITS=10 --PLIST=ParameterLists/ParameterList_GDSWCoarseOperator_Send.xml"
COMM mpi
NUM_MPI_PROCS 8
)
