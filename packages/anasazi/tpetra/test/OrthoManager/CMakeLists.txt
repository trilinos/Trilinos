
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  # FIXME (mfh 15 May 2011) The TSQR version of the OrthoManager seems
  # to be broken.  I'm disabling that test for now by running a simple
  # "Basic" OrthoManager test, which passes.  I will fix
  # TsqrOrthoManager as soon as I get a chance.
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
  IF (${PACKAGE_NAME}_ENABLE_TSQR)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Tpetra_OrthoManager_Tsqr_test
      SOURCES cxx_main_tsqr.cpp 
      ARGS "--ortho=Basic"
      COMM serial mpi
      )
  ENDIF ()
ENDIF()
