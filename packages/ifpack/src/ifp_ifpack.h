/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef _IFP_IFPACK_H_
#define _IFP_IFPACK_H_

#include "Ifpack_ConfigDefs.h"

#include <iostream>
#include <stdlib.h>
#include "ifp_arch.h"

#ifdef CRAY
extern "C" { void exit(int); }
#endif

#define ifp_error(msg,code) \
{ std::cerr << "IFPACK: " << msg << "    CODE " << code << std::endl; exit(1); }

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define ABS(x)   ((x)<0 ? (-(x)) : (x))
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#define SGN(x) ((x)<0.0 ? -1.0 : 1.0)

// miscellaneous prototypes for FORTRAN functions

extern "C"
{
void F77NAME(dgetri)(int *n, double *A, int *lda, int *ipiv,
    double *work, int *lwork, int *info);

void F77NAME(dgesvd) (char *, char *, int *, int *, const double *,
    int *, double *, double *, int *, double *, int *,
    double *, int *, int *);
}

#endif // _IFP_IFPACK_H_
