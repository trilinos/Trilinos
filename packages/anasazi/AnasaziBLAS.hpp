
#ifndef ANASAZI_BLAS_H
#define ANASAZI_BLAS_H

//! AnasaziBLAS:  The Anasazi BLAS Wrapper Class.
/*! The AnasaziBLAS class is a wrapper that encapsulates the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vector multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The AnasaziBLAS class provides C++ wrappers for the BLAS
    kernels in order to insulate the rest of Anasazi from the details of C++ to Fortran
    translation.
    A AnasaziBLAS object is essentially nothing, but allows access to the BLAS wrapper
    functions.
  
    AnasaziBLAS is a serial interface only.  This is appropriate since the standard 
    BLAS are only specified for serial execution (or shared memory parallel).
*/


class AnasaziBLAS {
    
  public:
  //@{ \name Constructors/Destructor.
   //! AnasaziBLAS Constructor.
  /*! Builds an instance of a serial BLAS object.
   */
  AnasaziBLAS(void);


  //! AnasaziBLAS Copy Constructor.
  /*! Makes an exact copy of an existing AnasaziBLAS instance.
  */
  AnasaziBLAS(const AnasaziBLAS& BLAS);

  //! AnasaziBLAS Destructor.
  virtual ~AnasaziBLAS(void);
  //@}
  
  //@{ \name Level 1 BLAS
  //! AnasaziBLAS one norm function (SASUM).
  float ASUM(int N, float * X) const;
  //! AnasaziBLAS one norm function (DASUM).
  double ASUM(int N, double * X) const;

  //! AnasaziBLAS dot product function (SDOT).
  float DOT(int N, float * X, float * Y) const;
  //! AnasaziBLAS dot product function (DDOT).
  double DOT(int N, double * X, double * Y) const;

  //! AnasaziBLAS norm function (SNRM2).
  float NRM2(int N, float * X) const;
  //! AnasaziBLAS norm function (DNRM2).
  double NRM2(int N, double * X) const;

  //! AnasaziBLAS vector scale function (SSCAL)
  void SCAL( int N, float ALPHA, float * X) const;
  //! AnasaziBLAS vector scale function (DSCAL)
  void SCAL( int N, double ALPHA, double * X) const;


  //! AnasaziBLAS arg maximum of absolute value function (ISAMAX)
  int IAMAX( int N, float * X) const;
  //! AnasaziBLAS arg maximum of absolute value function (IDAMAX)
  int IAMAX( int N, double * X) const;

  //! AnasaziBLAS vector update function (SAXPY)
  void AXPY( int N, float ALPHA, float * X, float * Y) const;
  //! AnasaziBLAS vector update function (DAXPY)
  void AXPY( int N, double ALPHA, double * X, double * Y) const;
  //@}

  //@{ \name Level 2 BLAS
  //! AnasaziBLAS matrix-vector multiply function (SGEMV)
  void GEMV(char TRANS, int M, int N,
         float ALPHA, float * A, int LDA, float * X,
         float BETA, float * Y) const;
  //! AnasaziBLAS matrix-vector multiply function (DGEMV)
  void GEMV(char TRANS, int M, int N,
         double ALPHA, double * A, int LDA, double * X,
         double BETA, double * Y) const;
  //@}


  //@{ \name Level 3 BLAS
  //! AnasaziBLAS matrix-matrix multiply function (SGEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const;
  //! AnasaziBLAS matrix-matrix multiply function (DGEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const;

  //! AnasaziBLAS symmetric matrix-matrix multiply function (SSYMM)
  void SYMM(char SIDE, char UPLO, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const;
  //! AnasaziBLAS matrix-matrix multiply function (DSYMM)
  void SYMM(char SIDE, char UPLO, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const;

  //! AnasaziBLAS triangular matrix-matrix multiply function (STRMM)
  void TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB) const;
  //! AnasaziBLAS triangular matrix-matrix multiply function (DTRMM)
  void TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB) const;
  //@}
};

// AnasaziBLAS constructor
inline AnasaziBLAS::AnasaziBLAS(void){}
// AnasaziBLAS constructor
inline AnasaziBLAS::AnasaziBLAS(const AnasaziBLAS& BLAS){}
// AnasaziBLAS destructor
inline AnasaziBLAS::~AnasaziBLAS(){}

#endif /* _ANASAZI_BLAS_H_ */
