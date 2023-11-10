TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IMEX_RK_Partitioned
  SOURCES Tempus_IMEX_RK_PartitionedTest.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE(
  IMEX_RK_Partitioned_Combined_FSA
  SOURCES Tempus_IMEX_RK_Partitioned_Combined_FSA.cpp Tempus_IMEX_RK_Partitioned_FSA.hpp
  TESTONLYLIBS tempus_test_models
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA
  NAME "IMEX_RK_Partitioned_Combined_FSA_Partitioned_IMEX_RK_1st_Order"
  ARGS "--method=\"Partitioned IMEX RK 1st order\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA
  NAME "IMEX_RK_Partitioned_Combined_FSA_Partitioned_IMEX_RK_SSP2"
  ARGS "--method=\"Partitioned IMEX RK SSP2\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA
  NAME "IMEX_RK_Partitioned_Combined_FSA_Partitioned_IMEX_RK_ARS_233"
  ARGS "--method=\"Partitioned IMEX RK ARS 233\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA
  NAME "IMEX_RK_Partitioned_Combined_FSA_General_Partioned_IMEX_RK"
  ARGS "--method=\"General Partitioned IMEX RK\""
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE(
  IMEX_RK_Partitioned_Combined_FSA_Tangent
  SOURCES Tempus_IMEX_RK_Partitioned_Combined_FSA_Tangent.cpp Tempus_IMEX_RK_Partitioned_FSA.hpp
  TESTONLYLIBS tempus_test_models
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Combined_FSA_Tangent_Partitioned_IMEX_RK_1st_Order"
  ARGS "--method=\"Partitioned IMEX RK 1st order\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Combined_FSA_Tangent_Partitioned_IMEX_RK_SSP2"
  ARGS "--method=\"Partitioned IMEX RK SSP2\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Combined_FSA_Tangent_Partitioned_IMEX_RK_ARS_233"
  ARGS "--method=\"Partitioned IMEX RK ARS 233\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Combined_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Combined_FSA_Tangent_General_Partioned_IMEX_RK"
  ARGS "--method=\"General Partitioned IMEX RK\""
  NUM_MPI_PROCS 1
  )



TRIBITS_ADD_EXECUTABLE(
  IMEX_RK_Partitioned_Staggered_FSA
  SOURCES Tempus_IMEX_RK_Partitioned_Staggered_FSA.cpp Tempus_IMEX_RK_Partitioned_FSA.hpp
  TESTONLYLIBS tempus_test_models
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Partitioned_IMEX_RK_1st_Order"
  ARGS "--method=\"Partitioned IMEX RK 1st order\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Partitioned_IMEX_RK_SSP2"
  ARGS "--method=\"Partitioned IMEX RK SSP2\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Partitioned_IMEX_RK_ARS_233"
  ARGS "--method=\"Partitioned IMEX RK ARS 233\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA
  NAME "IMEX_RK_Partitioned_Staggered_FSA_General_Partioned_IMEX_RK"
  ARGS "--method=\"General Partitioned IMEX RK\""
  NUM_MPI_PROCS 1
  )



TRIBITS_ADD_EXECUTABLE(
  IMEX_RK_Partitioned_Staggered_FSA_Tangent
  SOURCES Tempus_IMEX_RK_Partitioned_Staggered_FSA_Tangent.cpp Tempus_IMEX_RK_Partitioned_FSA.hpp
  TESTONLYLIBS tempus_test_models
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Tangent_Partitioned_IMEX_RK_1st_Order"
  ARGS "--method=\"Partitioned IMEX RK 1st order\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Tangent_Partitioned_IMEX_RK_SSP2"
  ARGS "--method=\"Partitioned IMEX RK SSP2\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Tangent_Partitioned_IMEX_RK_ARS_233"
  ARGS "--method=\"Partitioned IMEX RK ARS 233\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  IMEX_RK_Partitioned_Staggered_FSA_Tangent
  NAME "IMEX_RK_Partitioned_Staggered_FSA_Tangent_General_Partioned_IMEX_RK"
  ARGS "--method=\"General Partitioned IMEX RK\""
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_IMEX_RK_Partitioned_CopyFiles
  DEST_FILES Tempus_IMEX_RK_VanDerPol.xml
  EXEDEPS IMEX_RK_Partitioned IMEX_RK_Partitioned_Combined_FSA IMEX_RK_Partitioned_Combined_FSA_Tangent IMEX_RK_Partitioned_Staggered_FSA IMEX_RK_Partitioned_Staggered_FSA_Tangent
  )
