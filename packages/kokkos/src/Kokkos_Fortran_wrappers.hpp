//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_FORTRAN_WRAPPERS_H
#define KOKKOS_FORTRAN_WRAPPERS_H

#include "Kokkos_f77func.hpp"

#   define KOKKOS_DCRSMV_F77                 F77_FUNC_(kokkos_dcrsmv,KOKKOS_DCRSMV)

#ifdef __cplusplus
extern "C" {
#endif

void PREFIX KOKKOS_DCRSMV_F77(const int *itrans, const int *udiag, const int *numRows, 
			      const int *numCols, double values[],
			      int indices[], int profile[], double x[], double y[]);

#ifdef __cplusplus
}
#endif

#endif /* KOKKOS_FORTRAN_WRAPPERS_H */
