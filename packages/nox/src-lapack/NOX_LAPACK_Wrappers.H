// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef _NOX_LAPACK_WRAPPERS_H_
#define _NOX_LAPACK_WRAPPERS_H_

/* Declarations of LAPACK routines */

#ifdef __cplusplus
extern "C" {
#endif

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77 F77_FUNC(dgetrs,DGETRS)
#define DGEEV_F77 F77_FUNC(dgeev,DGEEV)

#ifdef HAVE_LAPACK_DGGEV
#define DGGEV_F77 F77_FUNC(dggev,DGGEV)
#endif
#ifdef HAVE_LAPACK_DGEGV
#define DGGEV_F77 F77_FUNC(dgegv,DGEGV)
#endif


/* B <- A \ B, A <- [L,U], piv <- P where P*A = L*U */
void DGESV_F77(const int *n, const int *nrhs, double *A, const int *lda, 
	       int *piv, double *B, const int *ldb, int *info);

/* A <- [L,U], piv <- P where P*A = L*U */
void DGETRF_F77(const int *m, const int *n, double *A, const int *lda, 
		int *piv, int *info);

/* B <- A \ B where A = [L,U], piv given from dgetrf_ */
void DGETRS_F77(const char *trans, const int *n, const int *nrhs, 
		const double *A, const int *lda, const int *piv,
		double *B, const int *ldb, int *info);

/* Compute left/right eigenvectors, eigenvalues */
void DGEEV_F77(const char *jobvl, const char *jobvr, const int *n,
	       double *A, const int *lda, double *wr, double *wi,
	       double *vl, const int *ldvl, double *vr, const int *ldvr,
	       double *work, const int *lwork, int *info);

#ifdef HAVE_LAPACK_GENEV
/* Compute left/right generalized eigenvectors, eigenvalues */
void DGGEV_F77(const char *jobvl, const char *jobvr, const int *n,
	       double *A, const int *lda, double *B, 
	       const int *ldb, double *alphar, double *alphai,
	       double *beta, double *vl, const int *ldvl, double *vr, 
	       const int *ldvr, double *work, const int *lwork, int *info);
#endif

#ifdef __cplusplus
}
#endif

#endif // _NOX_LAPACK_WRAPPERS_H_
