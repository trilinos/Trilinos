
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_LAPACK_H_
#define _EPETRA_LAPACK_H_

//! Epetra_LAPACK:  The Epetra LAPACK Wrapper Class.
/*! The Epetra_LAPACK class is a wrapper that encapsulates LAPACK
    (Linear Algebra Package).  LAPACK provides portable, high-
    performance implementations of linear, eigen, SVD, etc solvers.

    The standard LAPACK interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Epetra_LAPACK class provides C++ wrappers for the LAPACK
    kernels in order to insulate the rest of Epetra from the details of C++ to Fortran
    translation.
    A Epetra_LAPACK object is essentially nothing, but allows access to the LAPACK wrapper
    functions.
  
    Epetra_LAPACK is a serial interface only.  This is appropriate since the standard 
    LAPACK are only specified for serial execution (or shared memory parallel).
*/

#include "Epetra_Object.h"

class Epetra_LAPACK {
    
  public:
  //@{ \name Constructors/destructors.
  //! Epetra_LAPACK Constructor.
  /*! Builds an instance of a serial LAPACK object.
   */
  Epetra_LAPACK(void);


  //! Epetra_LAPACK Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_LAPACK instance.
  */
  Epetra_LAPACK(const Epetra_LAPACK& LAPACK);

  //! Epetra_LAPACK Destructor.
  virtual ~Epetra_LAPACK(void);
  //@}


  //@{ \name Symmetric Positive Definite linear system routines.
  
  //! Epetra_LAPACK factorization for positive definite matrix (SPOTRF)
  void POTRF( char UPLO, int N, float * A, int LDA, int * INFO) const;
  //! Epetra_LAPACK factorization for positive definite matrix (DPOTRF)
  void POTRF( char UPLO, int N, double * A, int LDA, int * INFO) const;

  //! Epetra_LAPACK solve (after factorization) for positive definite matrix (SPOTRS)
  void POTRS( char UPLO, int N, int NRHS, float * A, int LDA, float * X, int LDX, int * INFO) const;
  //! Epetra_LAPACK solve (after factorization) for positive definite matrix (DPOTRS)
  void POTRS( char UPLO, int N, int NRHS, double * A, int LDA, double * X, int LDX, int * INFO) const;

  //! Epetra_LAPACK inversion  for positive definite matrix (SPOTRI)
  void POTRI( char UPLO, int N, float * A, int LDA, int * INFO) const;
  //! Epetra_LAPACK inversion  for positive definite matrix (DPOTRI)
  void POTRI( char UPLO, int N, double * A, int LDA, int * INFO) const;

  //! Epetra_LAPACK condition number estimator for positive definite matrix (SPOCON)
  void POCON( char UPLO, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  //! Epetra_LAPACK condition number estimator for positive definite matrix (DPOCON)
  void POCON( char UPLO, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, int * INFO) const;

  //! Epetra_LAPACK factor and solve for positive definite matrix (SPOSV)
  void POSV( char UPLO, int N, int NRHS, float * A, int LDA, float * X, int LDX, int * INFO) const;
  //! Epetra_LAPACK factor and solve for positive definite matrix (DPOSV)
  void POSV( char UPLO, int N, int NRHS, double * A, int LDA, double * X, int LDX, int * INFO) const;

  //! Epetra_LAPACK equilibration for positive definite matrix (SPOEQU)
  void POEQU(int N, float * A, int LDA, float * S, float * SCOND, float * AMAX, int * INFO) const;
  //! Epetra_LAPACK equilibration for positive definite matrix (DPOEQU)
  void POEQU(int N, double * A, int LDA, double * S, double * SCOND, double * AMAX, int * INFO) const;

  //! Epetra_LAPACK solve driver for positive definite matrix (SPOSVX)
  void PORFS(char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     float * B, int LDB, float * X, int LDX, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Epetra_LAPACK solve driver for positive definite matrix (DPOSVX)
  void PORFS(char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! Epetra_LAPACK solve driver for positive definite matrix (SPOSVX)
  void POSVX(char FACT, char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     char EQUED, float * S, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Epetra_LAPACK solve driver for positive definite matrix (DPOSVX)
  void POSVX(char FACT, char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     char EQUED, double * S, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;
  //@}

  //@{ \name General linear system routines.

  //! Epetra_LAPACK simple driver to solve least-squares systems
  void GELS( char trans, int m, int n, int numrhs, double* a, int lda, 
	  double* b, int ldb, double* work, int lwork, int info) const;
  //! Epetra_LAPACK factorization for general matrix (SGETRF)
  void GETRF( int M, int N, float * A, int LDA, int * IPIV, int * INFO) const;
  //! Epetra_LAPACK factorization for general matrix (DGETRF)
  void GETRF( int M, int N, double * A, int LDA, int * IPIV, int * INFO) const;

  //! Epetra_LAPACK solve (after factorization) for general matrix (SGETRS)
  void GETRS( char TRANS, int N, int NRHS, float * A, int LDA, int * IPIV, float * X, int LDX, int * INFO) const;
  //! Epetra_LAPACK solve (after factorization) for general matrix (DGETRS)
  void GETRS( char TRANS, int N, int NRHS, double * A, int LDA, int * IPIV, double * X, int LDX, int * INFO) const;

  //! Epetra_LAPACK inversion  for general matrix (SGETRI)
  void GETRI( int N, float * A, int LDA, int * IPIV, float * WORK, int * LWORK, int * INFO) const;
  //! Epetra_LAPACK inversion  for general matrix (DGETRI)
  void GETRI( int N, double * A, int LDA, int * IPIV, double * WORK, int * LWORK, int * INFO) const;

  //! Epetra_LAPACK condition number estimator for general matrix (SGECON)
  void GECON( char NORM, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  //! Epetra_LAPACK condition number estimator for general matrix (DGECON)
  void GECON( char NORM, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, int * INFO) const;

  //! Epetra_LAPACK factor and solve for general matrix (SGESV)
  void GESV( int N, int NRHS, float * A, int LDA, int * IPIV, float * X, int LDX, int * INFO) const;
  //! Epetra_LAPACK factor and solve for general matrix (DGESV)
  void GESV( int N, int NRHS, double * A, int LDA, int * IPIV, double * X, int LDX, int * INFO) const;

  //! Epetra_LAPACK equilibration for general matrix (SGEEQU)
  void GEEQU(int M, int N, float * A, int LDA, float * R, float * C, float * ROWCND, float * COLCND, float * AMAX, int * INFO) const;
  //! Epetra_LAPACK equilibration for general matrix (DGEEQU)
  void GEEQU(int M, int N, double * A, int LDA, double * R, double * C, double * ROWCND, double * COLCND, double * AMAX, int * INFO) const;

  //! Epetra_LAPACK solve driver for general matrix (SGESVX)
  void GERFS(char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     int * IPIV, float * B, int LDB, float * X, int LDX, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Epetra_LAPACK solve driver for general matrix (DGESVX)
  void GERFS(char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     int * IPIV, double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;

  //! Epetra_LAPACK solve driver for general matrix (SGESVX)
  void GESVX(char FACT, char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, int * IPIV, 
	     char EQUED, float * R, float * C, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  //! Epetra_LAPACK solve driver for general matrix (DGESVX)
  void GESVX(char FACT, char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, int * IPIV, 
	     char EQUED, double * R, double * C, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;


  //! Epetra_LAPACK wrapper for reduction to Hessenberg form (SGEHRD)
  void GEHRD(int N, int ILO, int IHI, float * A, int LDA, float * TAU, float * WORK, int LWORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for reduction to Hessenberg form (DGEHRD)
  void GEHRD(int N, int ILO, int IHI, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Hessenberg routines
  //! Epetra_LAPACK wrapper for computing the eigenvalues of a real upper Hessenberg matrix (SHSEQR)
  void HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, float * H, int LDH, float * WR, float * WI,
	      float * Z, int LDZ, float * WORK, int LWORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for computing the eigenvalues of a real upper Hessenberg matrix (DHSEQR)
  void HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, double * H, int LDH, double * WR, double * WI,
	      double * Z, int LDZ, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Orthogonal matrix routines
  //! Epetra_LAPACK wrapper for generating a real orthogonal matrix Q defined by elementary reflectors. (SORGHR)
  void ORGHR( int N, int ILO, int IHI, float * A, int LDA, float * TAU, float * WORK, int LWORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for generating a real orthogonal matrix Q defined by elementary reflectors. (DORGHR)
  void ORGHR( int N, int ILO, int IHI, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) const;

  //! Epetra_LAPACK wrapper for applying an orthogonal matrix in-place (SORMHR)
  void ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, float * A, int LDA, 
	      float * TAU, float * C,
	      int LDC, float * WORK, int LWORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for applying an orthogonal matrix in-place (DORMHR)
  void ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, double * A, int LDA, 
	      double * TAU, double * C,
	      int LDC, double * WORK, int LWORK, int * INFO) const;
  //@}

  //@{ \name Triangular matrix routines

  //! Epetra_LAPACK wrapper for computing eigenvectors of a quasi-triangular/triagnular matrix (STREVC)
  /*! \warning HOWMNY = 'S" is not supported.
   */
  void TREVC( char SIDE, char HOWMNY, int * SELECT, int N, float * T, int LDT, float *VL, int LDVL,
	      float * VR, int LDVR, int MM, int * M, float * WORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for computing eigenvectors of a quasi-triangular/triagnular matrix (DTREVC)
  /*! \warning HOWMNY = 'S" is not supported.
   */
  void TREVC( char SIDE, char HOWMNY, int * SELECT, int N, double * T, int LDT, double *VL, int LDVL,
	      double * VR, int LDVR, int MM, int  *M, double * WORK, int * INFO) const;

  //! Epetra_LAPACK wrapper for reordering the real-Schur/Schur factorization of a matrix (STREXC)
  void TREXC( char COMPQ, int N, float * T, int LDT, float * Q, int LDQ, int IFST, int ILST, 
	      float * WORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for reordering the real-Schur/Schur factorization of a matrix (DTREXC)
  void TREXC( char COMPQ, int N, double * T, int LDT, double * Q, int LDQ, int IFST, int ILST, 
	      double * WORK, int * INFO) const;
  //@}

  //@{ \name Singular Value Decomposition matrix routines

  //! Epetra_LAPACK wrapper for computing the singular value decomposition (SGESVD)
  void GESVD( char JOBU, char JOBVT, int M, int N, float * A, int LDA, float * S, float * U,
	      int LDU, float * VT, int LDVT, float * WORK, int * LWORK, int * INFO) const;
  //! Epetra_LAPACK wrapper for computing the singular value decomposition (DGESVD)
  void GESVD( char JOBU, char JOBVT, int M, int N, double * A, int LDA, double * S, double * U,
	      int LDU, double * VT, int LDVT, double * WORK, int * LWORK, int * INFO) const;
   //@}
  //@{ \name Machine characteristics routines.
  //! Epetra_LAPACK wrapper for DLAMCH routine.  On out, T holds machine double precision floating point characteristics.  This information is returned by the Lapack routine.
  void LAMCH ( char CMACH, float & T) const;
  //! Epetra_LAPACK wrapper for SLAMCH routine.  On out, T holds machine single precision floating point characteristics.  This information is returned by the Lapack routine.
  void LAMCH ( char CMACH, double & T) const;
 //@}

};

// Epetra_LAPACK constructor
inline Epetra_LAPACK::Epetra_LAPACK(void){}
// Epetra_LAPACK constructor
inline Epetra_LAPACK::Epetra_LAPACK(const Epetra_LAPACK& LAPACK){}
// Epetra_LAPACK destructor
inline Epetra_LAPACK::~Epetra_LAPACK(){}

#endif /* _EPETRA_LAPACK_H_ */
