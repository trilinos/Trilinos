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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
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

#ifdef __cplusplus
}
#endif

#endif // _NOX_LAPACK_WRAPPERS_H_
