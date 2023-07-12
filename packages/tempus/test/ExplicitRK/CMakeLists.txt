TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK
  SOURCES Tempus_ExplicitRKTest.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE(
  ExplicitRK_Combined_FSA
  SOURCES Tempus_ExplicitRK_Combined_FSA.cpp Tempus_ExplicitRK_FSA.hpp
  TESTONLYLIBS tempus_test_models
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Forward_Euler"
  ARGS "--method=\"RK Forward Euler\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_RK_Explicit_4_Stage"
  ARGS "--method=\"RK Explicit 4 Stage\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_3_8_Rule"
  ARGS "--method=\"RK Explicit 3/8 Rule\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_4_Stage_3rd_Order_By_Runge"
  ARGS "--method=\"RK Explicit 4 Stage 3rd order by Runge\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_5_Stage_3rd_Order_By_Kinnmark_And_Gray"
  ARGS "--method=\"RK Explicit 5 Stage 3rd order by Kinnmark and Gray\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_3_Stage_3rd_Order"
  ARGS "--method=\"RK Explicit 3 Stage 3rd order\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_3_Stage_3rd_Order_TVD"
  ARGS "--method=\"RK Explicit 3 Stage 3rd order TVD\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_3_Stage_3rd_Order_By_Heun"
  ARGS "--method=\"RK Explicit 3 Stage 3rd order by Heun\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_Midpoint"
  ARGS "--method=\"RK Explicit Midpoint\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_RK_Explicit_Trapezoidal"
  ARGS "--method=\"RK Explicit Trapezoidal\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_Heuns_Method"
  ARGS "--method=\"Heuns Method\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Combined_FSA
  NAME "ExplicitRK_Combined_FSA_General_ERK"
  ARGS "--method=\"General ERK\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE(
  ExplicitRK_Staggered_FSA
  SOURCES Tempus_ExplicitRK_Staggered_FSA.cpp Tempus_ExplicitRK_FSA.hpp
  TESTONLYLIBS tempus_test_models
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Forward_Euler"
  ARGS "--method=\"RK Forward Euler\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_RK_Explicit_4_Stage"
  ARGS "--method=\"RK Explicit 4 Stage\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_3_8_Rule"
  ARGS "--method=\"RK Explicit 3/8 Rule\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_4_Stage_3rd_Order_By_Runge"
  ARGS "--method=\"RK Explicit 4 Stage 3rd order by Runge\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_5_Stage_3rd_Order_By_Kinnmark_And_Gray"
  ARGS "--method=\"RK Explicit 5 Stage 3rd order by Kinnmark and Gray\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_3_Stage_3rd_Order"
  ARGS "--method=\"RK Explicit 3 Stage 3rd order\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_3_Stage_3rd_Order_TVD"
  ARGS "--method=\"RK Explicit 3 Stage 3rd order TVD\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_3_Stage_3rd_Order_By_Heun"
  ARGS "--method=\"RK Explicit 3 Stage 3rd order by Heun\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_Midpoint"
  ARGS "--method=\"RK Explicit Midpoint\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_RK_Explicit_Trapezoidal"
  ARGS "--method=\"RK Explicit Trapezoidal\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_Heuns_Method"
  ARGS "--method=\"Heuns Method\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  ExplicitRK_Staggered_FSA
  NAME "ExplicitRK_Staggered_FSA_General_ERK"
  ARGS "--method=\"General ERK\""
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK_PseudoTransient_SA
  SOURCES Tempus_ExplicitRK_PseudoTransient_SA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK_ASA
  SOURCES Tempus_ExplicitRK_ASA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_ExplicitRK_CopyFiles
    DEST_FILES Tempus_ExplicitRK_SinCos.xml Tempus_ExplicitRK_SteadyQuadratic.xml Tempus_ExplicitRK_VanDerPol.xml
  EXEDEPS ExplicitRK
  )
