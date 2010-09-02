INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

#
# We do build MPI versions of the Combine and SequentialTsqr tests,
# but they just quiet all processes but Rank 0 and do the tests there.
#

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_Combine_Accuracy
  SOURCES Tsqr_TestCombine.cpp
  COMM serial mpi
  ARGS "--verify --nrows=1000 --ncols=5"
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_SequentialTsqr_Accuracy
  SOURCES Tsqr_TestSeqTsqr.cpp
  COMM serial mpi
  ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_DistTsqr_Accuracy
  SOURCES Tsqr_TestDistTsqr.cpp
  COMM serial mpi
  ARGS "--verify --ncols=10"
  STANDARD_PASS_OUTPUT
)

IF (TSQR_ENABLE_Intel_TBB)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TSQR_TbbTsqr_Accuracy
    SOURCES Tsqr_TestTbbTsqr.cpp
    COMM serial mpi
    ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   TSQR_MpiSeqTsqr_Accuracy
#   SOURCES Tsqr_TestMpiSeqTsqr.cpp
#   COMM mpi
#   ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
#   STANDARD_PASS_OUTPUT
# )

