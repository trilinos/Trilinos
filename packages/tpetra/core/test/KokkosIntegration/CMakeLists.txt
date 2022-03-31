
ASSERT_DEFINED (Tpetra_ENABLE_CUDA)

IF (Tpetra_ENABLE_CUDA AND BUILD_SHARED_LIBS)
  MESSAGE(STATUS "Tpetra: Enabling KokkosIntegration Tests")

  TRIBITS_ADD_LIBRARY(
    Intercept
    SOURCES Intercept.cpp
    SHARED
    TESTONLY
  )

  TRIBITS_ADD_LIBRARY(
    ApiTest
    SOURCES ApiTest.cpp
    SHARED
    TESTONLY
  )

  TRIBITS_ADD_EXECUTABLE(
    KokkosUVM
    SOURCES KokkosUVM.cpp 
    TESTONLYLIBS ApiTest
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    KokkosDeepcopy
    SOURCES KokkosDeepcopy.cpp
    TESTONLYLIBS ApiTest
    COMM serial mpi
  )

 TRIBITS_ADD_TEST(
   KokkosDeepcopy
   NAME KokkosDeepcopy
   COMM mpi
   NUM_MPI_PROCS 1
   ENVIRONMENT LD_PRELOAD=./libIntercept.so
   STANDARD_PASS_OUTPUT
 )


ELSE()
MESSAGE(STATUS "Tpetra: Disabling KokkosIntegration Tests")
ENDIF()
