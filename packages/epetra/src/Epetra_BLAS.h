
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRA_BLAS_H
#define EPETRA_BLAS_H

#include "Epetra_Object.h"

//! Epetra_BLAS:  The Epetra BLAS Wrapper Class.
/*! The Epetra_BLAS class is a wrapper that encapsulates the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vectoer multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Epetra_BLAS class provides C++ wrappers for the BLAS
    kernels in order to insulate the rest of Epetra from the details of C++ to Fortran
    translation.
    A Epetra_BLAS object is essentially nothing, but allows access to the BLAS wrapper
    functions.
  
    Epetra_BLAS is a serial interface only.  This is appropriate since the standard 
    BLAS are only specified for serial execution (or shared memory parallel).
*/


class Epetra_BLAS {
    
  public:
  //@{ \name Constructors/Destructor.
   //! Epetra_BLAS Constructor.
  /*! Builds an instance of a serial BLAS object.
   */
  Epetra_BLAS(void);


  //! Epetra_BLAS Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_BLAS instance.
  */
  Epetra_BLAS(const Epetra_BLAS& BLAS);

  //! Epetra_BLAS Destructor.
  virtual ~Epetra_BLAS(void);
  //@}
  
  //@{ \name Level 1 BLAS
  //! Epetra_BLAS one norm function (SASUM).
  float ASUM(int N, float * X, int INCX = 1) const;
  //! Epetra_BLAS one norm function (DASUM).
  double ASUM(int N, double * X, int INCX = 1) const;

  //! Epetra_BLAS dot product function (SDOT).
  float DOT(int N, float * X, float * Y, int INCX = 1, int INCY = 1) const;
  //! Epetra_BLAS dot product function (DDOT).
  double DOT(int N, double * X, double * Y, int INCX = 1, int INCY = 1) const;

  //! Epetra_BLAS norm function (SNRM2).
  float NRM2(int N, float * X, int INCX = 1) const;
  //! Epetra_BLAS norm function (DNRM2).
  double NRM2(int N, double * X, int INCX = 1) const;

  //! Epetra_BLAS vector scale function (SSCAL)
  void SCAL( int N, float ALPHA, float * X, int INCX = 1) const;
  //! Epetra_BLAS vector scale function (DSCAL)
  void SCAL( int N, double ALPHA, double * X, int INCX = 1) const;

  //! Epetra_BLAS vector copy function (SCOPY)
  void COPY( int N, float * X, float * Y, int INCX = 1, int INCY = 1) const;
  //! Epetra_BLAS vector scale function (DCOPY)
  void COPY( int N, double * X, double * Y, int INCX = 1, int INCY = 1) const;

  //! Epetra_BLAS arg maximum of absolute value function (ISAMAX)
  int IAMAX( int N, float * X, int INCX = 1) const;
  //! Epetra_BLAS arg maximum of absolute value function (IDAMAX)
  int IAMAX( int N, double * X, int INCX = 1) const;

  //! Epetra_BLAS vector update function (SAXPY)
  void AXPY( int N, float ALPHA, float * X, float * Y, int INCX = 1, int INCY = 1) const;
  //! Epetra_BLAS vector update function (DAXPY)
  void AXPY( int N, double ALPHA, double * X, double * Y, int INCX = 1, int INCY = 1) const;
  //@}

  //@{ \name Level 2 BLAS
  //! Epetra_BLAS matrix-vector multiply function (SGEMV)
  void GEMV(char TRANS, int M, int N,
         float ALPHA, float * A, int LDA, float * X,
         float BETA, float * Y, int INCX = 1, int INCY = 1) const;
  //! Epetra_BLAS matrix-vector multiply function (DGEMV)
  void GEMV(char TRANS, int M, int N,
         double ALPHA, double * A, int LDA, double * X,
         double BETA, double * Y, int INCX = 1, int INCY = 1) const;
  //@}


  //@{ \name Level 3 BLAS
  //! Epetra_BLAS matrix-matrix multiply function (SGEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const;
  //! Epetra_BLAS matrix-matrix multiply function (DGEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const;

  //! Epetra_BLAS symmetric matrix-matrix multiply function (SSYMM)
  void SYMM(char SIDE, char UPLO, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const;
  //! Epetra_BLAS matrix-matrix multiply function (DSYMM)
  void SYMM(char SIDE, char UPLO, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const;

  //! Epetra_BLAS triangular matrix-matrix multiply function (STRMM)
  void TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB) const;
  //! Epetra_BLAS triangular matrix-matrix multiply function (DTRMM)
  void TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB) const;
  //@}
};

// Epetra_BLAS constructor
inline Epetra_BLAS::Epetra_BLAS(void){}
// Epetra_BLAS constructor
inline Epetra_BLAS::Epetra_BLAS(const Epetra_BLAS& BLAS){}
// Epetra_BLAS destructor
inline Epetra_BLAS::~Epetra_BLAS(){}

#endif /* EPETRA_BLAS_H */
