
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_EXECUTABLE(
    krino_unit
    SOURCES ${SOURCES}
    COMM serial mpi
    NOEXEPREFIX NOEXESUFFIX)
    
TRIBITS_ADD_TEST(krino_unit
    NAME Krino_Basic
    ARGS "--gtest_filter=-*fuzz*:*Random*:HamiltonJacobi.*:DecomposeMesh*Sensitivities*:*Refinement.performance*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)

TRIBITS_ADD_TEST(krino_unit
    NAME Krino_Fuzz
    ARGS "--gtest_filter=*fuzz*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
    
TRIBITS_ADD_TEST(krino_unit
    NAME Krino_Robustness
    ARGS "--gtest_filter=*Random*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
    
TRIBITS_ADD_TEST(krino_unit
    NAME Krino_HamiltonJacobi
    ARGS "--gtest_filter=HamiltonJacobi.*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
    
TRIBITS_ADD_TEST(krino_unit
    NAME Krino_Topology_Optimization
    ARGS "--gtest_filter=DecomposeMesh*Sensitivities*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
    
TRIBITS_ADD_TEST(krino_unit
    NAME Krino_Refinement_Performance
    ARGS "--gtest_filter=*Refinement.performance*"
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX NOEXESUFFIX)
