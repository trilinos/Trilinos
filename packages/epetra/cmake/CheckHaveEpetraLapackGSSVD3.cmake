# LAPACK 3.6.0 deprecates dggsvd and sggsvd (Issue #480)
# The new versions of dggsvd and sggsvd are called dggsvd3 and sggsvd3
#
# dggsvd (JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO)
# sggsvd (JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO)
#
# dggsvd3 (JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO)
# sggsvd3 (JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO)
#
# The new parameter is LWORK
#
#   Information is at:
#   http://www.netlib.org/lapack/explore-html/d6/db3/dggsvd3_8f_a4a187519e5c71da3b3f67c85e9baf0f2.html#a4a187519e5c71da3b3f67c85e9baf0f2

INCLUDE(CheckFunctionExists)

FUNCTION(CHECK_HAVE_EPETRA_LAPACK_GSSVD3 VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_LAPACK_LIBRARIES} ${TPL_BLAS_LIBRARIES})

  CHECK_FUNCTION_EXISTS( dggsvd3   HAVE_dggsvd3 )
  CHECK_FUNCTION_EXISTS( dggsvd3_  HAVE_dggsvd3_POST )
  CHECK_FUNCTION_EXISTS( DGGSVD3   HAVE_DGGSVD3 )
  CHECK_FUNCTION_EXISTS( DGGSVD3_  HAVE_DGGSVD3_POST )

  IF( HAVE_dggsvd3 OR HAVE_dggsvd3_POST OR HAVE_DGGSVD3 OR HAVE_DGGSVD3_POST )
    MESSAGE( "Found new version of lapack. dggsvd3 is available." )

    # Note calling CHECK_FUNCTION_EXISTS( dggsvd3_, ${VARNAME} ) is ok for scope
    # but setting like this after several checks, it seems PARENT_SCOPE must be added
    SET( ${VARNAME} 1 PARENT_SCOPE )

  ELSE()
    MESSAGE( "Did not find new version of lapack. dggsvd3 is not available." )
  ENDIF()

ENDFUNCTION()
