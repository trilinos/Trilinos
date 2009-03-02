
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ml_Validate
  SOURCES cxx_main.cpp
  COMM serial mpi
  )

PACKAGE_ADD_TEST(
  ml_Validate
  ARGS
    "--num-time-steps=10 --soln-verb-level=VERB_EXTREME"
    "--params-file=simpleAdjoint_n4.xml --num-time-steps=10
--soln-verb-level=VERB_EXTREME"
  STANDARD_PASS_OUTPUT
  COMM serial mpi  
  NUM_MPI_PROCS 1 
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(ValidateParamsCopyFiles
  DEST_FILES input1.xml baseline1.txt evaluateCriteria
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/ValidateParams
  EXEDEPS ml_Validate
  )
