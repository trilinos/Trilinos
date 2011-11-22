
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MLP_Sys
  SOURCES MultiLevelPreconditioner_Sys.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(MLP_Sys_Copy
  SOURCE_FILES c1.dat c2.dat c3.dat elasticity.dat
  EXEDEPS MLP_Sys
  )
