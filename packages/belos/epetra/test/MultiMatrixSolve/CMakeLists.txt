
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bl_gmres_multidiag
  SOURCES test_bl_gmres_multidiag.cpp BelosLinearMultiShiftProblem.hpp 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pseudo_gmres_multi_hb
    SOURCES test_pseudo_gmres_multi_hb.cpp BelosLinearMultiShiftProblem.hpp
    COMM serial mpi
    ARGS "--verbose --filename=orsirr1.hb"
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestGmresMultiFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES orsirr1.hb
    EXEDEPS pseudo_gmres_multi_hb
  )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)    
