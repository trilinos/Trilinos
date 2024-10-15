
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef AMESOS_SCALAPACK_WRAPPERS_H
#define AMESOS_SCALAPACK_WRAPPERS_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Epetra_ConfigDefs.h"
#include "Epetra_LAPACK_wrappers.h"
#include <stdio.h>
#include <string.h>

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#if defined(CRAY_T3X)

#define SL_INIT_F77  F77_FUNC_(sl_init,SL_INIT)
#define BLACS_GRIDINFO_F77  F77_FUNC_(blacs_gridinfo,BLACS_GRIDINFO)
#define PDGETRF_F77  F77_FUNC(psgetrf,PSGETRF)
#define PDGETRS_F77  F77_FUNC(psgetrs,PSGETRS)
#define DESCINIT_F77  F77_FUNC(descinit,DESCINIT)

#endif
#if defined(INTEL_CXML)

#define SL_INIT_F77  F77_FUNC_(sl_init,SL_INIT)
#define BLACS_GRIDINFO_F77  F77_FUNC_(blacs_gridinfo,BLACS_GRIDINFO)
#define PDGETRF_F77  F77_FUNC(pdgetrf,PDGETRF)
#define PDGETRS_F77  F77_FUNC(pdgetrs,PDGETRS)
#define DESCINIT_F77  F77_FUNC(descinit,DESCINIT)

#endif
#if defined(INTEL_MKL)

#define SL_INIT_F77  F77_FUNC_(sl_init,SL_INIT_)
#define BLACS_GRIDINFO_F77  F77_FUNC_(blacs_gridinfo,BLACS_GRIDINFO_)
#define PDGETRF_F77  F77_FUNC(pdgetrf,PDGETRF)
#define PDGETRS_F77  F77_FUNC(pdgetrs,PDGETRS)
#define DESCINIT_F77  F77_FUNC(descinit,DESCINIT)

#endif

#else

/* Use autoconf's definition of F77_FUNC 
   unless using old make system */

#define SL_INIT_F77  F77_FUNC_(sl_init,SL_INIT)
#define BLACS_GRIDINFO_F77  F77_FUNC_(blacs_gridinfo,BLACS_GRIDINFO)

#define PDGETRF_F77  F77_FUNC(pdgetrf,PDGETRF)
#define PDGETRS_F77  F77_FUNC(pdgetrs,PDGETRS)
#define DESCINIT_F77  F77_FUNC(descinit,DESCINIT)

#endif

#ifdef __cplusplus
extern "C" {
#endif

  /* ScaLAPACK and BLACS initialization routines */
  void PREFIX SL_INIT_F77(int* blacs_context, const int* nprow, const int* npcol);
  void PREFIX DESCINIT_F77(int *DescA, const int* m, const int* n, const int* mblock, 
			   const int* nblock, const int* rsrc, const int* csrc, const int* blacs_context,
			   const int* Lda, int* ierr);
  void PREFIX BLACS_GRIDINFO_F77(int* blacs_context, const int* nprow, const int* npcol,
				 const int* myrow, const int* mycol);
  /* Double precision ScaLAPACK linear solvers */
  void PREFIX PDGETRF_F77(const int* m, const int* n, double* A, const int* Ai, const int* Aj, 
			  const int* DescA, int* ipiv, int* info);
  void PREFIX PDGETRS_F77(Epetra_fcd, const int* n, const int* nrhs, 
			  const double* A, const int* Ai, const int* Aj, 
			  const int* DescA, const int* ipiv, double* X, const int* Xi, const int* Xj,
			  const int* DescX, int* info);

#ifdef __cplusplus
}
#endif

#endif /* AMESOS_SCALAPACK_WRAPPERS_H */
