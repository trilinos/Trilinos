
#ifndef ANASAZI_LAPACK_H
#define ANASAZI_LAPACK_H

//! AnasaziLAPACK:  The Anasazi LAPACK Wrapper Class.
/*! The AnasaziLAPACK class is a wrapper that encapsulates LAPACK
    (Linear Algebra Package).  LAPACK provides portable, high-
    performance implementations of linear, eigen, SVD, etc solvers.

    The standard LAPACK interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The AnasaziLAPACK class provides C++ wrappers for the LAPACK
    kernels in order to insulate the rest of Anasazi from the details of C++ to Fortran
    translation.
    A AnasaziLAPACK object is essentially nothing, but allows access to the LAPACK wrapper
    functions.
  
    AnasaziLAPACK is a serial interface only.  This is appropriate since the standard 
    LAPACK are only specified for serial execution (or shared memory parallel).

    This wrapper class originated from Epetra, and was adapted for Anasazi.
*/

class AnasaziLAPACK {
    
  public:
  //@{ \name Constructors/destructors.
  //! AnasaziLAPACK Constructor.
  /*! Builds an instance of a serial LAPACK object.
   */
  AnasaziLAPACK(void);


  //! AnasaziLAPACK Copy Constructor.
  /*! Makes an exact copy of an existing AnasaziLAPACK instance.
  */
  AnasaziLAPACK(const AnasaziLAPACK& LAPACK);

  //! AnasaziLAPACK Destructor.
  virtual ~AnasaziLAPACK(void);
  //@}


  //@{ \name Symmetric Positive Definite linear system routines.
  
  //! AnasaziLAPACK factorization for positive definite matrix (SPOTRF)
  void POTRF( char UPLO, int N, float * A, int LDA, int * INFO) const;
  //! AnasaziLAPACK factorization for positive definite matrix (DPOTRF)
  void POTRF( char UPLO, int N, double * A, int LDA, int * INFO) const;

  //! AnasaziLAPACK solve (after factorization) for positive definite matrix (SPOTRS)
  void POTRS( char UPLO, int N, int NRHS, float * A, int LDA, float * X, int LDX, int * INFO) const;
  //! AnasaziLAPACK solve (after factorization) for positive definite matrix (DPOTRS)
  void POTRS( char UPLO, int N, int NRHS, double * A, int LDA, double * X, int LDX, int * INFO) const;

  //! AnasaziLAPACK inversion  for positive definite matrix (SPOTRI)
  void POTRI( char UPLO, int N, float * A, int LDA, int * INFO) const;
  //! AnasaziLAPACK inversion  for positive definite matrix (DPOTRI)
  void POTRI( char UPLO, int N, double * A, int LDA, int * INFO) const;

  //! AnasaziLAPACK condition number estimator for positive definite matrix (SPOCON)
  void POCON( char UPLO, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  //! AnasaziLAPACK condition number estimator for positive definite matrix (DPOCON)
  void POCON( char UPLO, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, int * INFO) const;

  //! AnasaziLAPACK factor and solve for positive definite matrix (SPOSV)
  void POSV( char UPLO, int N, int NRHS, float * A, int LDA, float * X, int LDX, int * INFO) const;
  //! AnasaziLAPACK factor and solve for positive definite matrix (DPOSV)
  void POSV( char UPLO, int N, int NRHS, double * A, int LDA, double * X, int LDX, int * INFO) const;

  //! AnasaziLAPACK equilibration for positive definite matrix (SPOEQU)
  void POEQU(int N, float * A, int LDA, float * S, float * SCOND, float * AMAX, int * INFO) const;
  //! AnasaziLAPACK equilibration for positive definite matrix (DPOEQU)
  void POEQU(int N, double * A, int LDA, double * S, double * SCOND, double * AMAX, int * INFO) const;

  //! AnasaziLAPACK solve driver for positive definite matrix (SPOSVX)
  void PORFS(char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     float * B, int LDB, float * X, int LDX, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! AnasaziLAPACK solve driver for positive definite matrix (DPOSVX)
  void PORFS(char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! AnasaziLAPACK solve driver for positive definite matrix (SPOSVX)
  void POSVX(char FACT, char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     char EQUED, float * S, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! AnasaziLAPACK solve driver for positive definite matrix (DPOSVX)
  void POSVX(char FACT, char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     char EQUED, double * S, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;
  //@}

  //@{ \name General linear system routines.

  //! AnasaziLAPACK simple driver to solve least-squares systems
  void GELS( char trans, int m, int n, int numrhs, double* a, int lda, 
	  double* b, int ldb, double* work, int lwork, int info) const;
  //! AnasaziLAPACK factorization for general matrix (SGETRF)
  void GETRF( int M, int N, float * A, int LDA, int * IPIV, int * INFO) const;
  //! AnasaziLAPACK factorization for general matrix (DGETRF)
  void GETRF( int M, int N, double * A, int LDA, int * IPIV, int * INFO) const;

  //! AnasaziLAPACK solve (after factorization) for general matrix (SGETRS)
  void GETRS( char TRANS, int N, int NRHS, float * A, int LDA, int * IPIV, float * X, int LDX, int * INFO) const;
  //! AnasaziLAPACK solve (after factorization) for general matrix (DGETRS)
  void GETRS( char TRANS, int N, int NRHS, double * A, int LDA, int * IPIV, double * X, int LDX, int * INFO) const;

  //! AnasaziLAPACK inversion  for general matrix (SGETRI)
  void GETRI( int N, float * A, int LDA, int * IPIV, float * WORK, int * LWORK, int * INFO) const;
  //! AnasaziLAPACK inversion  for general matrix (DGETRI)
  void GETRI( int N, double * A, int LDA, int * IPIV, double * WORK, int * LWORK, int * INFO) const;

  //! AnasaziLAPACK condition number estimator for general matrix (SGECON)
  void GECON( char NORM, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  //! AnasaziLAPACK condition number estimator for general matrix (DGECON)
  void GECON( char NORM, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, int * INFO) const;

  //! AnasaziLAPACK factor and solve for general matrix (SGESV)
  void GESV( int N, int NRHS, float * A, int LDA, int * IPIV, float * X, int LDX, int * INFO) const;
  //! AnasaziLAPACK factor and solve for general matrix (DGESV)
  void GESV( int N, int NRHS, double * A, int LDA, int * IPIV, double * X, int LDX, int * INFO) const;

  //! AnasaziLAPACK equilibration for general matrix (SGEEQU)
  void GEEQU(int M, int N, float * A, int LDA, float * R, float * C, float * ROWCND, float * COLCND, float * AMAX, int * INFO) const;
  //! AnasaziLAPACK equilibration for general matrix (DGEEQU)
  void GEEQU(int M, int N, double * A, int LDA, double * R, double * C, double * ROWCND, double * COLCND, double * AMAX, int * INFO) const;

  //! AnasaziLAPACK solve driver for general matrix (SGESVX)
  void GERFS(char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     int * IPIV, float * B, int LDB, float * X, int LDX, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! AnasaziLAPACK solve driver for general matrix (DGESVX)
  void GERFS(char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     int * IPIV, double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! AnasaziLAPACK solve driver for general matrix (SGESVX)
  void GESVX(char FACT, char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, int * IPIV, 
	     char EQUED, float * R, float * C, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! AnasaziLAPACK solve driver for general matrix (DGESVX)
  void GESVX(char FACT, char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, int * IPIV, 
	     char EQUED, double * R, double * C, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! AnasaziLAPACK eigenproblem driver (Schur form) for nonsymmetric matrix (SGEES)
  void GEES(char JOBVS, char SORT, int * SELECT, int N, float * A, int LDA, int SDIM, 
	     float* WR, float * WI, float* VS, int LDVS, float* WORK, int LWORK, 
	     int *BWORK, int* INFO ) const; 
  //! AnasaziLAPACK eigenproblem driver (Schur form) for nonsymmetric matrix (DGEES)
  void GEES(char JOBVS, char SORT, int * SELECT, int N, double * A, int LDA, int SDIM, 
	     double* WR, double * WI, double* VS, int LDVS, double* WORK, int LWORK, 
	     int *BWORK, int* INFO ) const; 

  //! AnasaziLAPACK wrapper for reduction to Hessenberg form (SGEHRD)
  void GEHRD(int N, int ILO, int IHI, float * A, int LDA, float * TAU, float * WORK, int LWORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for reduction to Hessenberg form (DGEHRD)
  void GEHRD(int N, int ILO, int IHI, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Hessenberg routines
  //! AnasaziLAPACK wrapper for computing the eigenvalues of a real upper Hessenberg matrix (SHSEQR)
  void HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, float * H, int LDH, float * WR, float * WI,
	      float * Z, int LDZ, float * WORK, int LWORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for computing the eigenvalues of a real upper Hessenberg matrix (DHSEQR)
  void HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, double * H, int LDH, double * WR, double * WI,
	      double * Z, int LDZ, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Orthogonal matrix routines
  //! AnasaziLAPACK wrapper for generating a real orthogonal matrix Q defined by elementary reflectors. (SORGHR)
  void ORGHR( int N, int ILO, int IHI, float * A, int LDA, float * TAU, float * WORK, int LWORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for generating a real orthogonal matrix Q defined by elementary reflectors. (DORGHR)
  void ORGHR( int N, int ILO, int IHI, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) const;

  //! AnasaziLAPACK wrapper for applying an orthogonal matrix in-place (SORMHR)
  void ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, float * A, int LDA, 
	      float * TAU, float * C,
	      int LDC, float * WORK, int LWORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for applying an orthogonal matrix in-place (DORMHR)
  void ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, double * A, int LDA, 
	      double * TAU, double * C,
	      int LDC, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Triangular matrix routines

  //! AnasaziLAPACK wrapper for computing eigenvectors of a quasi-triangular matrix (STREVC)
  /*! \warning HOWMNY = 'S" is not supported.
   */
  void TREVC( char SIDE, char HOWMNY, int * SELECT, int N, float * T, int LDT, float *VL, int LDVL,
	      float * VR, int LDVR, int MM, int * M, float * WORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for computing eigenvectors of a quasi-triangular matrix (DTREVC)
  /*! \warning HOWMNY = 'S" is not supported.
   */
  void TREVC( char SIDE, char HOWMNY, int * SELECT, int N, double * T, int LDT, double *VL, int LDVL,
	      double * VR, int LDVR, int MM, int  *M, double * WORK, int * INFO) const;

  //! AnasaziLAPACK wrapper for reordering the real-Schur/Schur factorization of a matrix (STREXC)
  void TREXC( char COMPQ, int N, float * T, int LDT, float * Q, int LDQ, int IFST, int ILST, 
	      float * WORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for reordering the real-Schur/Schur factorization of a matrix (DTREXC)
  void TREXC( char COMPQ, int N, double * T, int LDT, double * Q, int LDQ, int IFST, int ILST, 
	      double * WORK, int * INFO) const;
  //@}

  //@{ \name Singular Value Decomposition matrix routines

  //! AnasaziLAPACK wrapper for computing the singular value decomposition (SGESVD)
  void GESVD( char JOBU, char JOBVT, int M, int N, float * A, int LDA, float * S, float * U,
	      int LDU, float * VT, int LDVT, float * WORK, int * LWORK, int * INFO) const;
  //! AnasaziLAPACK wrapper for computing the singular value decomposition (DGESVD)
  void GESVD( char JOBU, char JOBVT, int M, int N, double * A, int LDA, double * S, double * U,
	      int LDU, double * VT, int LDVT, double * WORK, int * LWORK, int * INFO) const;
   //@}

  //@{ \name Machine characteristics routines.
  //! AnasaziLAPACK wrapper for DLAMCH routine that returns machine double precision floating point characteristics
  double DLAMCH( char CMACH) const; 
  //! AnasaziLAPACK wrapper for SLAMCH routine that returns machine single precision floating point characteristics
  float SLAMCH( char CMACH) const;
 //@}

  //@{ \name Auxiliary routines
  //! AnasaziLAPACK wrapper for DLAPY2 routine that returns sqrt(x**2 + y**2), taking care not to cause unnecessary overflow.
  double LAPY2( double X, double Y ) const;
  //! AnasaziLAPACK wrapper for SLAPY2 routine that returns sqrt(x**2 + y**2), taking care not to cause unnecessary overflow.
  float LAPY2( float X, float Y ) const;  
  //@}
};

// AnasaziLAPACK constructor
inline AnasaziLAPACK::AnasaziLAPACK(void){}
// AnasaziLAPACK constructor
inline AnasaziLAPACK::AnasaziLAPACK(const AnasaziLAPACK& LAPACK){}
// AnasaziLAPACK destructor
inline AnasaziLAPACK::~AnasaziLAPACK(){}

#endif /* ANASAZI_LAPACK_H */
