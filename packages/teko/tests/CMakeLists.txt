INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
  testdriver
  SOURCES
      ./src/Epetra/tEpetraOperatorWrapper.cpp
      ./src/Epetra/tEpetraOperatorWrapper.hpp
      ./src/Epetra/tEpetraThyraConverter.cpp
      ./src/Epetra/tEpetraThyraConverter.hpp
      ./src/Epetra/tInterlacedEpetra.cpp
      ./src/Epetra/tInterlacedEpetra.hpp
      ./src/Epetra/tStridedEpetraOperator.cpp
      ./src/Epetra/tStridedEpetraOperator.hpp
      ./src/tBlockJacobiPreconditionerFactory.cpp
      ./src/tBlockJacobiPreconditionerFactory.hpp
      ./src/tBlockLowerTriInverseOp.cpp
      ./src/tBlockLowerTriInverseOp.hpp
      ./src/tBlockUpperTriInverseOp.cpp
      ./src/tBlockUpperTriInverseOp.hpp
      ./src/tExplicitOps.cpp
      ./src/tExplicitOps.hpp
      ./src/tGraphLaplacian.cpp
      ./src/tGraphLaplacian.hpp
      ./src/tJacobi2x2PreconditionerFactory.cpp
      ./src/tJacobi2x2PreconditionerFactory.hpp
      ./src/tLSCHIntegrationTest.cpp
      ./src/tLSCHIntegrationTest.hpp
      ./src/tLSCIntegrationTest.cpp
      ./src/tLSCIntegrationTest.hpp
      ./src/tLSCStabilized.cpp
      ./src/tLSCStabilized.hpp
      ./src/tLSCStablePreconditionerFactory.cpp
      ./src/tLSCStablePreconditionerFactory.hpp
      ./src/tLU2x2PreconditionerFactory.cpp
      ./src/tLU2x2PreconditionerFactory.hpp
      ./src/tParallelInverse.cpp
      ./src/tParallelInverse.hpp
      ./src/tSIMPLEPreconditionerFactory.cpp
      ./src/tSIMPLEPreconditionerFactory.hpp
      ./Test_Utils.cpp
      ./Test_Utils.hpp
      ./testdriver.cpp
  COMM mpi serial
  )

PACKAGE_ADD_TEST(
   testdriver
   ARGS --fast
   NUM_MPI_PROCS 4
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

PACKAGE_ADD_TEST(
   testdriver
   ARGS --fast
   NUM_MPI_PROCS 1
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

PACKAGE_COPY_FILES_TO_BINARY_DIR(testdriver_copyfiles
  SOURCE_FILES lsc_B_2.mm
     lsc_Bt_2.mm
     lsc_F_2.mm
     lsc_Qu_2.mm
     lsc_exact_2.mm
     lsc_rhs.mm
     nsjac.mm
     nsjac_1.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
  EXEDEPS testdriver
  )
