INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
IF(Teuchos_ENABLE_COMPLEX)
  TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

  #
  # mfh 20 Jan 2014: This test ensures the following:
  #
  # 1. Belos::GCRODRSolMgr can compile whether its ScalarType (first)
  #    template parameter is real or complex.
  # 2. Belos::GCRODRSolMgr's constructor throws std::logic_error if and
  #    only if its ScalarType (first) template parameter is complex.
  #
  # At some point, if GCRODRSolMgr gets fixed so that it works with
  # complex ScalarType, the second test will no longer pass.  This will
  # be a good thing!  The test should still be built in that case, in
  # order to demonstrate that GCRODRSolMgr compiles for complex
  # ScalarType.  However, in that case, the TEST_THROW macro should be
  # changed to TEST_NOTHROW, and the macro's second argument should be
  # removed.
  # 
  # This test requires that Trilinos was compiled with complex
  # arithmetic support enabled.
  #
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_gcrodr_complex
    SOURCES
      test_gcrodr_complex.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    )

  #
  # mfh 20 Jan 2014: This test ensures the following:
  #
  # 1. Belos::RCGSolMgr can compile whether its ScalarType (first)
  #    template parameter is real or complex.
  # 2. Belos::RCGSolMgr's constructor throws std::logic_error if and
  #    only if its ScalarType (first) template parameter is complex.
  #
  # At some point, if RCGSolMgr gets fixed so that it works with
  # complex ScalarType, the second test will no longer pass.  This will
  # be a good thing!  The test should still be built in that case, in
  # order to demonstrate that RCGSolMgr compiles for complex
  # ScalarType.  However, in that case, the TEST_THROW macro should be
  # changed to TEST_NOTHROW, and the macro's second argument should be
  # removed.
  # 
  # This test requires that Trilinos was compiled with complex
  # arithmetic support enabled.
  #
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_rcg_complex
    SOURCES
      test_rcg_complex.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    )

  #
  # mfh 20 Jan 2014: This test ensures the following:
  #
  # 1. Belos::LSQRSolMgr can compile whether its ScalarType (first)
  #    template parameter is real or complex.
  # 2. Belos::LSQRSolMgr's constructor throws std::logic_error if and
  #    only if its ScalarType (first) template parameter is complex.
  #
  # At some point, if LSQRSolMgr gets fixed so that it works with
  # complex ScalarType, the second test will no longer pass.  This will
  # be a good thing!  The test should still be built in that case, in
  # order to demonstrate that LSQRSolMgr compiles for complex
  # ScalarType.  However, in that case, the TEST_THROW macro should be
  # changed to TEST_NOTHROW, and the macro's second argument should be
  # removed.
  # 
  # This test requires that Trilinos was compiled with complex
  # arithmetic support enabled.
  #
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_lsqr_complex
    SOURCES
      test_lsqr_complex.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    )

  #
  # mfh 20 Jan 2014: This test ensures that Belos::SolverFactory can
  # compile whether its ScalarType (first) template parameter is real
  # or complex.
  #
  # This test requires that Trilinos was compiled with complex
  # arithmetic support enabled.
  #
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_factory_complex
    SOURCES
      test_factory_complex.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    )

ENDIF(Teuchos_ENABLE_COMPLEX)
