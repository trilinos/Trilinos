#ifndef _PETRA_BLAS_H_
#define _PETRA_BLAS_H_

//! Petra_BLAS:  The Petra BLAS Wrapper Class.
/*! The Petra_BLAS class is a wrapper that encapsulates the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vectoer multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Petra_BLAS class provides C++ wrappers for the BLAS
    kernels in order to insulate the rest of Petra from the details of C++ to Fortran
    translation.
    A Petra_BLAS object is essentially nothing, but allows access to the BLAS wrapper
    functions.
  
    Petra_BLAS is a serial interface only.  This is appropriate since the standard 
    BLAS are only specified for serial execution (or shared memory parallel).
*/

#include "Petra_Petra.h"

class Petra_BLAS {
    
  public:
  //! Petra_BLAS Constructor.
  /*! Builds an instance of a serial BLAS object.
   */
  Petra_BLAS(void);


  //! Petra_BLAS Copy Constructor.
  /*! Makes an exact copy of an existing Petra_BLAS instance.
  */
  Petra_BLAS(const Petra_BLAS& BLAS);

  //! Petra_BLAS Destructor.
  virtual ~Petra_BLAS(void);
  
  //! Petra_BLAS one norm function (SASUM).
  float ASUM(int N, float * X) const;
  //! Petra_BLAS one norm function (DASUM).
  double ASUM(int N, double * X) const;

  //! Petra_BLAS dot product function (SDOT).
  float DOT(int N, float * X, float * Y) const;
  //! Petra_BLAS dot product function (DDOT).
  double DOT(int N, double * X, double * Y) const;

  //! Petra_BLAS norm function (SNRM2).
  float NRM2(int N, float * X) const;
  //! Petra_BLAS norm function (DNRM2).
  double NRM2(int N, double * X) const;

  //! Petra_BLAS vector scale function (SSCAL)
  void SCAL( int N, float ALPHA, float * X) const;
  //! Petra_BLAS vector scale function (DSCAL)
  void SCAL( int N, double ALPHA, double * X) const;


  //! Petra_BLAS arg maximum of absolute value function (ISAMAX)
  int IAMAX( int N, float * X) const;
  //! Petra_BLAS arg maximum of absolute value function (IDAMAX)
  int IAMAX( int N, double * X) const;

  //! Petra_BLAS vector update function (SAXPY)
  void AXPY( int N, float ALPHA, float * X, float * Y) const;
  //! Petra_BLAS vector update function (DAXPY)
  void AXPY( int N, double ALPHA, double * X, double * Y) const;

  //! Petra_BLAS matrix-vector multiply function (SGEMV)
  void GEMV(char TRANS, int M, int N,
         float ALPHA, float * A, int LDA, float * X,
         float BETA, float * Y) const;
  //! Petra_BLAS matrix-vector multiply function (DGEMV)
  void GEMV(char TRANS, int M, int N,
         double ALPHA, double * A, int LDA, double * X,
         double BETA, double * Y) const;

  //! Petra_BLAS matrix-matrix multiply function (SGEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const;
  //! Petra_BLAS matrix-matrix multiply function (DGEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const;

  //! Petra_BLAS triangular matrix-matrix multiply function (STRMM)
  void TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB) const;
  //! Petra_BLAS triangular matrix-matrix multiply function (DTRMM)
  void TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB) const;
};

// Petra_BLAS constructor
inline Petra_BLAS::Petra_BLAS(void){}
// Petra_BLAS constructor
inline Petra_BLAS::Petra_BLAS(const Petra_BLAS& BLAS){}
// Petra_BLAS destructor
inline Petra_BLAS::~Petra_BLAS(){}

#endif /* _PETRA_BLAS_H_ */
