INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Experimental)
IF(${PACKAGE_NAME}_ENABLE_Experimental)
  # 
  # This test depends on experimental code.
  #
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    projected_least_squares_solver
    SOURCES test_belos_projected_least_squares_solver.cpp
    COMM serial mpi
    )
ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental
