
# As of <= 14 Mar 2012, MultiVectorFiller doesn't work yet with
# multiple MPI processes.  Thus, I'll leave it out of the MPI tests
# for now.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVectorFiller_SerialNodeTest
  SOURCES 
    MultiVectorFiller_SerialNodeTest
  COMM serial
  STANDARD_PASS_OUTPUT
  )

IF (TpetraCore_ENABLE_TpetraClassic AND TpetraClassic_ENABLE_TBB AND ${PACKAGE_NAME}_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_TbbTest
  #   SOURCES MultiVectorFiller_TbbTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   )
ENDIF()

IF (TpetraCore_ENABLE_TpetraClassic AND TpetraClassic_ENABLE_ThreadPool AND ${PACKAGE_NAME}_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_ThreadPoolTest
  #   SOURCES MultiVectorFiller_ThreadPoolTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   )
ENDIF()

IF (TpetraCore_ENABLE_TpetraClassic AND TpetraClassic_ENABLE_Thrust AND TpetraClassic_ENABLE_CUDA_DOUBLE AND ${PACKAGE_NAME}_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_ThrustTest
  #   SOURCES MultiVectorFiller_ThrustTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   )
ENDIF()

