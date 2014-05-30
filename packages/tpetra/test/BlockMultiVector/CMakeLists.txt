
INCLUDE(TribitsAddExecutableAndTest)

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_COMPILE(BLOCKMV_UNIT_TEST_OBJECTS_CUDA BlockMultiVector_UnitTests_CU.cu)
  CUDA_COMPILE(EXPBLOCKMULTIVECTOR_OBJECTS_CUDA ExpBlockMultiVector_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockMultiVector_UnitTests
  SOURCES
    BlockMultiVector_UnitTests
    ${BLOCKMV_UNIT_TEST_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockMultiVector
  SOURCES
    ExpBlockMultiVector
    ${EXPBLOCKMULTIVECTOR_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

# mfh 23 May 2014: Temporarily disabled test due to test (not class!)
# build error with Scalar=std::complex.  Please fix test and reenable
# ASAP.
#
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   ExpBlockCrsMatrix
#   SOURCES
#     ExpBlockCrsMatrix
#     ${EXPBLOCKMULTIVECTOR_OBJECTS_CUDA}
#     ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#   )
