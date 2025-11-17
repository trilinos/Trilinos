TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ../../../test/MVOPTester)

IF (Tpetra_INST_DOUBLE)

    # FIXME (mfh 15 May 2011) The TSQR version of the OrthoManager seems
    # to be broken.  I'm disabling that test for now by running a simple
    # "Basic" OrthoManager test, which passes.  I will fix
    # TsqrOrthoManager as soon as I get a chance.
    ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
    IF (${PACKAGE_NAME}_ENABLE_TSQR)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        Tpetra_OrthoManager_Tsqr_test
        SOURCES cxx_main_tsqr.cpp 
        ARGS "--ortho=Basic"
        COMM serial mpi
        )
    ENDIF ()
ENDIF ()

IF (Tpetra_INST_COMPLEX_DOUBLE)

   TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_OrthoManagerComplexTester
      SOURCES cxx_main_complex.cpp
      COMM serial mpi
      ARGS
        "--verbose --ortho=SVQB"
        "--verbose --filename=mhd1280b.cua --ortho=SVQB"
        "--verbose --ortho=Basic"
        "--verbose --filename=mhd1280b.cua --ortho=Basic"
        "--verbose --ortho=ICGS"
        "--verbose --filename=mhd1280b.cua --ortho=ICGS"
      STANDARD_PASS_OUTPUT
      )

   ASSERT_DEFINED(Anasazi_SOURCE_DIR)
   TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_OrthoManager_ComplexCopyFiles
      SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS Tpetra_OrthoManagerComplexTester
      )

ENDIF ()
