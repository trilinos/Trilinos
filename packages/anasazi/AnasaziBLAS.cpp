
/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
*/

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, 1
#else
#define CHAR_MACRO(char_var) &char_var
#endif

#include "AnasaziBLAS.hpp"
#include "AnasaziBLASwrappers.hpp"


//=============================================================================
float AnasaziBLAS::ASUM(int N, float * X) const {
  int one = 1;
 return(SASUM_F77(&N, X, &one));
}
//=============================================================================
double AnasaziBLAS::ASUM(int N, double * X) const {
  int one = 1;
  return(DASUM_F77(&N, X, &one));
}
//=============================================================================
float AnasaziBLAS::DOT(int N, float * X, float * Y) const {
  int one = 1;
  return(SDOT_F77(&N, X, &one, Y, &one));
}
//=============================================================================
double AnasaziBLAS::DOT(int N, double * X, double * Y) const {
  int one = 1;
  return(DDOT_F77(&N, X, &one, Y, &one));
}
//=============================================================================
float AnasaziBLAS::NRM2(int N, float * X) const {
  int one = 1;
  return(SNRM2_F77(&N, X, &one));
}
//=============================================================================
double AnasaziBLAS::NRM2(int N, double * X) const {
  int one = 1;
  return(DNRM2_F77(&N, X, &one));
}
//=============================================================================
void AnasaziBLAS::SCAL(int N, float ALPHA, float * X) const {
  int one = 1;
  SSCAL_F77(&N, &ALPHA, X, &one);
  return;
}
//=============================================================================
void AnasaziBLAS::SCAL(int N, double ALPHA, double * X) const {
  int one = 1;
  DSCAL_F77(&N, &ALPHA, X, &one);
  return;
}
//=============================================================================
int AnasaziBLAS::IAMAX(int N, float * X) const {
  int one = 1;
  return(ISAMAX_F77(&N, X, &one)-1);// Note that we return base zero result
}
//=============================================================================
int AnasaziBLAS::IAMAX(int N, double * X) const {
  int one = 1;
  return(IDAMAX_F77(&N, X, &one)-1);// Note that we return base zero result
}
//=============================================================================
void AnasaziBLAS::AXPY(int N, float ALPHA, float * X, float * Y) const {
  int one = 1;
  SAXPY_F77(&N, &ALPHA, X, &one, Y, &one);
}
//=============================================================================
void AnasaziBLAS::AXPY(int N, double ALPHA, double * X, double * Y) const {
  int one = 1;
  DAXPY_F77(&N, &ALPHA, X, &one, Y, &one);
}
//=============================================================================
void AnasaziBLAS::GEMV(char TRANS, int M, int N,
		      float ALPHA, float * A, int LDA, float * X,
		      float BETA, float * Y) const {
  int one = 1;
  SGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
}
//=============================================================================
void AnasaziBLAS::GEMV(char TRANS, int M, int N,
		      double ALPHA, double * A, int LDA, double * X,
		      double BETA, double * Y) const {
  int one = 1;
  DGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
}

//=============================================================================
void AnasaziBLAS::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const {

  SGEMM_F77(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}

//=============================================================================
void AnasaziBLAS::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const {

  DGEMM_F77(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, &ALPHA,
	 A, &LDA, B, &LDB, &BETA, C, &LDC);
}
//=============================================================================
void AnasaziBLAS::SYMM(char SIDE, char UPLO, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const {

  SSYMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), &M, &N, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}

//=============================================================================
void AnasaziBLAS::SYMM(char SIDE, char UPLO, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const {

  DSYMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), &M, &N, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}
//=============================================================================
void AnasaziBLAS::TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB) const {

  STRMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANSA), CHAR_MACRO(DIAG), 
	  &M, &N, &ALPHA, A, &LDA, B, &LDB);
}
//=============================================================================
void AnasaziBLAS::TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB) const {

  DTRMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANSA), CHAR_MACRO(DIAG), 
	    &M, &N, &ALPHA, A, &LDA, B, &LDB);
}
