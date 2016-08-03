# LAPACK 3.6.0 deprecates dggsvd and sggsvd (Issue #480)
# The new versions of dggsvd and sggsvd are called dggsvd3 and sggsvd3
# This probes for the deprecated version
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

INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_HAVE_LAPACK_GSSVD3 VARNAME)

  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_LAPACK_LIBRARIES})

  SET(SOURCE
  "

#define F77_BLAS_MANGLE${F77_BLAS_MANGLE}

#define DGGSVD_F77   F77_BLAS_MANGLE(dggsvd3,DGGSVD3)

extern \"C\" { void DGGSVD_F77(
        const char * JOBU,
        const char * JOBV,
        const char * JOBQ,
        const int * M,
        const int * N,
        const int * P,
        int * K,
        int * L,
        double * A,
        const int * LDA,
        double * B,
        const int * LDB,
        double * ALPHA,
        double * BETA,
        double * U,
        const int * LDU,
        double * V,
        const int * LDV,
        double * Q,
        const int * LDQ,
        double * WORK,
        const int * LWORK, // the new parameter LWORK
        int * IWORK,
        int * INFO); }

int main()
{
  const char JOBU = 'N';
  const char JOBV = 'N';
  const char JOBQ = 'N';
  const int M = 2;
  const int N = 2;
  const int P = 2;
  int K = 0;
  int L = 0;
  double A[4] = { 1.0, 1.0, 1.0, 1.0 };
  int LDA = 2;
  double B[4] = { 1.0, 1.0, 1.0, 1.0 };
  int LDB = 2;
  double alpha[2] = { 1.0, 1.0 };
  double beta[2] = { 1.0, 1.0 };
  double U[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int LDU = 2;
  double V[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int LDV = 2;
  double Q[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int LDQ = 2;
  double WORK[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int LWORK = 4;
  int IWORK[2] = { 1, 1 };
  int INFO = -9999; // will check for valid exit code 0

  DGGSVD_F77(&JOBU, &JOBV, &JOBQ,
    &M, &N, &P,
    &K, &L,
    A, &LDA, B, &LDB,
    alpha, beta,
    U, &LDU, V, &LDV, Q, &LDQ,
    WORK,
    &LWORK,   // the new parameter LWORK
    IWORK, &INFO);

  if( INFO == 0 )
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
  "
  )

  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})

ENDFUNCTION()
