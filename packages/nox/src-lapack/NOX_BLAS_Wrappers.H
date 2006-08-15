// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef _NOX_BLAS_WRAPPERS_H_
#define _NOX_BLAS_WRAPPERS_H_

/* Declerations of Legacy BLAS routines */

#ifdef __cplusplus
extern "C" {
#endif

#define DAXPY_F77 F77_FUNC(daxpy,DAXPY)
#define DDOT_F77 F77_FUNC(ddot,DDOT)
#define DSCAL_F77 F77_FUNC(dscal,DSCAL)
#define DNRM2_F77 F77_FUNC(dnrm2,DNRM2)
#define DASUM_F77 F77_FUNC(dasum,DASUM)
#define DGEMV_F77 F77_FUNC(dgemv,DGEMV)
#define DGEMM_F77 F77_FUNC(dgemm,DGEMM)

/* y <- a*x + y */
void DAXPY_F77(const int *n, const double *a, const double *x, 
	       const int *incx, double *y, const int *incy);

/* v <- x'*y */
double DDOT_F77(const int *n, const double *x, const int *incx, 
		const double *y, const int *incy);

/* x <- a*x */
void DSCAL_F77(const int *n, const double *a, double *x, 
	       const int *incx);

/*  v <- ||x||_2 */
double DNRM2_F77(const int *n, const double *x, const int *incx);

/* v <- ||x||_1 */
double DASUM_F77(const int *n, const double *x, const int *incx);

/* y <- a*op(A)*x + b*y where op(A) = A or op(A) = A' */
void DGEMV_F77(const char *T, const int *m, const int *n, const double *a,
	       const double *A, const int *lda, const double *x, 
	       const int *incx, const double *b, double *y, const int *incy);

/* C <- a*op(A)*op(B) + b*C  where op(X) = X or op(X) = X' */
void DGEMM_F77(const char *Ta, const char *Tb, const int *m, const int *n, 
	       const int *k, const double *a, const double *A, const int *lda, 
	       const double *B, const int *ldb, const double *b, 
	       double *C, const int *ldc);

#ifdef __cplusplus
}
#endif

#endif //_NOX_BLAS_WRAPPERS_H_
