

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
IF(Teuchos_ENABLE_COMPLEX)

  TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    bl_cg_complex_hb
    SOURCES test_bl_cg_complex_hb.cpp
    ARGS
    "--verbose --filename=mhd1280b.cua"
    "--verbose --filename=mhd1280b.cua --pseudo"
    "--verbose --filename=mhd1280b.cua --use-single-red"
    XHOST KALLIKRATES
  )

  IF(NOT ${${PROJECT_NAME}_HOSTNAME} STREQUAL "KALLIKRATES")
    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGComplexFiles
      SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS bl_cg_complex_hb
    )
  ENDIF()

ENDIF(Teuchos_ENABLE_COMPLEX)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bl_cg_real_hb
  SOURCES test_bl_cg_real_hb.cpp
  ARGS
  "--verbose --filename=A.hb"
  "--verbose --filename=A.hb --pseudo"
  "--verbose --filename=A.hb --use-single-red"
  "--verbose --filename=A.hb --combine-conv-inner"
  "--verbose --filename=A.hb --pseudo --combine-conv-inner"
  "--verbose --filename=A.hb --use-single-red --combine-conv-inner"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGFiles2
  SOURCE_FILES A.hb
)
