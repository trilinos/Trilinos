// 27-May-2002 General cleanup. Changed method names to fit namingConvention (nothing changed).

#ifndef _TPETRA_LAPACK_H_
#define _TPETRA_LAPACK_H_

//! Tpetra::LAPACK:  The Tpetra LAPACK Wrapper Class.
/*! The Tpetra::LAPACK class is a wrapper that encapsulates LAPACK
    (Linear Algebra Package).  LAPACK provides portable, high-
    performance implementations of linear, eigen, SVD, etc solvers.

    The standard LAPACK interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Tpetra::LAPACK class provides C++ wrappers for the LAPACK
    kernels in order to insulate the rest of Tpetra from the details of C++
    to Fortran translation.

    A Tpetra::LAPACK object is essentially nothing, but allows access to the LAPACK wrapper
    functions. It should also be noted that at present, Tpetra::LAPACK is an almost exact
    duplicate of Epetra_LAPACK.
  
    Tpetra::LAPACK is a serial interface only.  This is appropriate since the standard 
    LAPACK are only specified for serial execution (or shared memory parallel).
*/

#include "Tpetra_Object.h"


namespace Tpetra
{
class LAPACK
{    
  public:
  //@{ \name Constructors/destructors.
  //! Tpetra::LAPACK Constructor.
  /*! Builds an instance of a serial LAPACK object.
   */
  inline LAPACK(void){};


  //! Tpetra::LAPACK Copy Constructor.
  /*! Makes an exact copy of an existing Tpetra::LAPACK instance.
  */
  inline LAPACK(const LAPACK& LAPACK){};

  //! Tpetra::LAPACK Destructor.
  inline virtual ~LAPACK(void){};
  //@}


  //@{ \name Symmetric Positive Definite linear system routines.
  
  //! Tpetra::LAPACK factorization for positive definite matrix (SPOTRF)
  void POTRF( char UPLO, int N, float * A, int LDA, int * INFO) const;
  //! Tpetra::LAPACK factorization for positive definite matrix (DPOTRF)
  void POTRF( char UPLO, int N, double * A, int LDA, int * INFO) const;

  //! Tpetra::LAPACK solve (after factorization) for positive definite matrix (SPOTRS)
  void POTRS( char UPLO, int N, int NRHS, float * A, int LDA, float * X, int LDX, int * INFO) const;
  //! Tpetra::LAPACK solve (after factorization) for positive definite matrix (DPOTRS)
  void POTRS( char UPLO, int N, int NRHS, double * A, int LDA, double * X, int LDX, int * INFO) const;

  //! Tpetra::LAPACK inversion  for positive definite matrix (SPOTRI)
  void POTRI( char UPLO, int N, float * A, int LDA, int * INFO) const;
  //! Tpetra::LAPACK inversion  for positive definite matrix (DPOTRI)
  void POTRI( char UPLO, int N, double * A, int LDA, int * INFO) const;

  //! Tpetra::LAPACK condition number estimator for positive definite matrix (SPOCON)
  void POCON( char UPLO, int N, float * A, int LDA, float ANORM, float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  //! Tpetra::LAPACK condition number estimator for positive definite matrix (DPOCON)
  void POCON( char UPLO, int N, double * A, int LDA, double ANORM, double * RCOND, double * WORK, int * IWORK, int * INFO) const;

  //! Tpetra::LAPACK factor and solve for positive definite matrix (SPOSV)
  void POSV( char UPLO, int N, int NRHS, float * A, int LDA, float * X, int LDX, int * INFO) const;
  //! Tpetra::LAPACK factor and solve for positive definite matrix (DPOSV)
  void POSV( char UPLO, int N, int NRHS, double * A, int LDA, double * X, int LDX, int * INFO) const;

  //! Tpetra::LAPACK equilibration for positive definite matrix (SPOEQU)
  void POEQU(int N, float * A, int LDA, float * S, float * SCOND, float * AMAX, int * INFO) const;
  //! Tpetra::LAPACK equilibration for positive definite matrix (DPOEQU)
  void POEQU(int N, double * A, int LDA, double * S, double * SCOND, double * AMAX, int * INFO) const;

  //! Tpetra::LAPACK solve driver for positive definite matrix (SPOSVX)
  void PORFS(char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, float * B, int LDB, float * X, 
	     int LDX, float * FERR,float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Tpetra::LAPACK solve driver for positive definite matrix (DPOSVX)
  void PORFS(char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, double * B, int LDB, double * X, 
	     int LDX, double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! Tpetra::LAPACK solve driver for positive definite matrix (SPOSVX)
  void POSVX(char FACT, char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, char EQUED, float * S, float * B, 
	     int LDB, float * X, int LDX, float * RCOND, float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Tpetra::LAPACK solve driver for positive definite matrix (DPOSVX)
  void POSVX(char FACT, char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, char EQUED, double * S, double * B, 
	     int LDB, double * X, int LDX, double * RCOND, double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;
  //@}

  //@{ \name General linear system routines.

  //! Tpetra::LAPACK simple driver to solve least-squares systems
  void GELS( char trans, int m, int n, int numrhs, double* a, int lda, 
	  double* b, int ldb, double* work, int lwork, int info) const;
  //! Tpetra::LAPACK factorization for general matrix (SGETRF)
  void GETRF( int M, int N, float * A, int LDA, int * IPIV, int * INFO) const;
  //! Tpetra::LAPACK factorization for general matrix (DGETRF)
  void GETRF( int M, int N, double * A, int LDA, int * IPIV, int * INFO) const;

  //! Tpetra::LAPACK solve (after factorization) for general matrix (SGETRS)
  void GETRS( char TRANS, int N, int NRHS, float * A, int LDA, int * IPIV, float * X, int LDX, int * INFO) const;
  //! Tpetra::LAPACK solve (after factorization) for general matrix (DGETRS)
  void GETRS( char TRANS, int N, int NRHS, double * A, int LDA, int * IPIV, double * X, int LDX, int * INFO) const;

  //! Tpetra::LAPACK inversion  for general matrix (SGETRI)
  void GETRI( int N, float * A, int LDA, int * IPIV, float * WORK, int * LWORK, int * INFO) const;
  //! Tpetra::LAPACK inversion  for general matrix (DGETRI)
  void GETRI( int N, double * A, int LDA, int * IPIV, double * WORK, int * LWORK, int * INFO) const;

  //! Tpetra::LAPACK condition number estimator for general matrix (SGECON)
  void GECON( char NORM, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  //! Tpetra::LAPACK condition number estimator for general matrix (DGECON)
  void GECON( char NORM, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, int * INFO) const;

  //! Tpetra::LAPACK factor and solve for general matrix (SGESV)
  void GESV( int N, int NRHS, float * A, int LDA, int * IPIV, float * X, int LDX, int * INFO) const;
  //! Tpetra::LAPACK factor and solve for general matrix (DGESV)
  void GESV( int N, int NRHS, double * A, int LDA, int * IPIV, double * X, int LDX, int * INFO) const;

  //! Tpetra::LAPACK equilibration for general matrix (SGEEQU)
  void GEEQU(int M, int N, float * A, int LDA, float * R, float * C, float * ROWCND, float * COLCND, float * AMAX, int * INFO) const;
  //! Tpetra::LAPACK equilibration for general matrix (DGEEQU)
  void GEEQU(int M, int N, double * A, int LDA, double * R, double * C, double * ROWCND, double * COLCND, double * AMAX, int * INFO) const;

  //! Tpetra::LAPACK solve driver for general matrix (SGESVX)
  void GERFS(char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     int * IPIV, float * B, int LDB, float * X, int LDX, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Tpetra::LAPACK solve driver for general matrix (DGESVX)
  void GERFS(char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     int * IPIV, double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! Tpetra::LAPACK solve driver for general matrix (SGESVX)
  void GESVX(char FACT, char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, int * IPIV, 
	     char EQUED, float * R, float * C, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Tpetra::LAPACK solve driver for general matrix (DGESVX)
  void GESVX(char FACT, char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, int * IPIV, 
	     char EQUED, double * R, double * C, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;


  //! Tpetra::LAPACK wrapper for reduction to Hessenberg form (SGEHRD)
  void GEHRD(int N, int ILO, int IHI, float * A, int LDA, float * TAU, float * WORK, int LWORK, int * INFO) const;
  //! Tpetra::LAPACK wrapper for reduction to Hessenberg form (DGEHRD)
  void GEHRD(int N, int ILO, int IHI, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Hessenberg routines
  //! Tpetra::LAPACK wrapper for computing the eigenvalues of a real upper Hessenberg matrix (SHSEQR)
  void HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, float * H, int LDH, float * WR, float * WI,
	      float * Z, int LDZ, float * WORK, int LWORK, int * INFO) const;
  //! Tpetra::LAPACK wrapper for computing the eigenvalues of a real upper Hessenberg matrix (DHSEQR)
  void HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, double * H, int LDH, double * WR, double * WI,
	      double * Z, int LDZ, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Orthogonal matrix routines
  //! Tpetra::LAPACK wrapper for generating a real orthogonal matrix Q defined by elementary reflectors. (SORGHR)
  void ORGHR( int N, int ILO, int IHI, float * A, int LDA, float * TAU, float * WORK, int LWORK, int * INFO) const;
  //! Tpetra::LAPACK wrapper for generating a real orthogonal matrix Q defined by elementary reflectors. (DORGHR)
  void ORGHR( int N, int ILO, int IHI, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) const;

  //! Tpetra::LAPACK wrapper for applying an orthogonal matrix in-place (SORMHR)
  void ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, float * A, int LDA, 
	      float * TAU, float * C,
	      int LDC, float * WORK, int LWORK, int * INFO) const;
  //! Tpetra::LAPACK wrapper for applying an orthogonal matrix in-place (DORMHR)
  void ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, double * A, int LDA, 
	      double * TAU, double * C,
	      int LDC, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Triangular matrix routines

  //! Tpetra::LAPACK wrapper for computing eigenvectors of a quasi-triangular/triagnular matrix (STREVC)
  /*! \warning HOWMNY = 'S" is not supported.
   */
  void TREVC( char SIDE, char HOWMNY, int * SELECT, int N, float * T, int LDT, float *VL, int LDVL,
	      float * VR, int LDVR, int MM, int * M, float * WORK, int * INFO) const;
  //! Tpetra::LAPACK wrapper for computing eigenvectors of a quasi-triangular/triagnular matrix (DTREVC)
  /*! \warning HOWMNY = 'S" is not supported.
   */
  void TREVC( char SIDE, char HOWMNY, int * SELECT, int N, double * T, int LDT, double *VL, int LDVL,
	      double * VR, int LDVR, int MM, int  *M, double * WORK, int * INFO) const;

  //! Tpetra::LAPACK wrapper for reordering the real-Schur/Schur factorization of a matrix (STREXC)
  void TREXC( char COMPQ, int N, float * T, int LDT, float * Q, int LDQ, int IFST, int ILST, 
	      float * WORK, int * INFO) const;
  //! Tpetra::LAPACK wrapper for reordering the real-Schur/Schur factorization of a matrix (DTREXC)
  void TREXC( char COMPQ, int N, double * T, int LDT, double * Q, int LDQ, int IFST, int ILST, 
	      double * WORK, int * INFO) const;
  //@}

  //@{ \name Machine characteristics routines.
  //! Tpetra::LAPACK wrapper for DLAMCH routine that returns machine double precision floating point characteristics
  double DLAMCH( char CMACH) const; 
  //! Tpetra::LAPACK wrapper for SLAMCH routine that returns machine single precision floating point characteristics
  float SLAMCH( char CMACH) const;
 //@}

};
#include "Tpetra_LAPACK.cpp"
} // namespace Tpetra
#endif /* _TPETRA_LAPACK_H_ */
