INCLUDE(TribitsAddExecutableAndTest)   
TRIBITS_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultMultiply  
  SOURCES 
    CrsMatrix_DefaultMultiply.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(  CrsMatrix_DefaultSolve  
  SOURCES 
    CrsMatrix_DefaultSolve.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 

# not the most important thing at the moment, will come back to it
# TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_UnitTests 
#   SOURCES   
#     CrsMatrix_UnitTests.cpp 
#     ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
#   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
# ) 

#IF (Kokkos_ENABLE_ThreadPool) 
#  TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_FirstTouchTPI
#    SOURCES   
#      CrsMatrix_FirstTouchTPI.cpp
#      ${TEUCHOS_STD_UNIT_TEST_MAIN}
#    COMM serial mpi
#    NUM_MPI_PROCS 1
#    STANDARD_PASS_OUTPUT
#    DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#    INSTALLABLE
#  ) 
#ENDIF (Kokkos_ENABLE_ThreadPool) 

#IF (Kokkos_ENABLE_TBB) 
#  TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_FirstTouchTBB
#    SOURCES   
#      CrsMatrix_FirstTouchTBB.cpp
#      ${TEUCHOS_STD_UNIT_TEST_MAIN}
#    COMM serial mpi
#    NUM_MPI_PROCS 1
#    STANDARD_PASS_OUTPUT
#    DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#    INSTALLABLE
#  ) 
#ENDIF (Kokkos_ENABLE_TBB) 

#IF (Kokkos_ENABLE_CUSP) 
#  TRIBITS_ADD_EXECUTABLE_AND_TEST( CrsMatrix_CUSP
#    SOURCES   
#      CrsMatrix_CUSP.cpp
#      ${TEUCHOS_STD_UNIT_TEST_MAIN}
#    COMM serial mpi
#    NUM_MPI_PROCS 1
#    STANDARD_PASS_OUTPUT
#    DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#    INSTALLABLE
#  ) 
#ENDIF()
