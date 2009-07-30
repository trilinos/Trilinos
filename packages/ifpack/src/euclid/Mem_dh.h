/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef MEM_DH_DH
#define MEM_DH_DH

#include "euclid_common.h"
#ifdef __cplusplus
extern "C" {
#endif

extern void  Mem_dhCreate(Mem_dh *m);
extern void  Mem_dhDestroy(Mem_dh m);

extern void *Mem_dhMalloc(Mem_dh m, size_t size);
extern void  Mem_dhFree(Mem_dh m, void *ptr);
  /* preceeding two are called via the macros
     MALLOC_DH and FREE_DH; see "euclid_config.h"
   */

extern void  Mem_dhPrint(Mem_dh m, FILE* fp, bool allPrint);
  /* prints memory usage statistics;  "allPrint" is only
     meaningful when running in MPI mode.
   */

#ifdef __cplusplus
}
#endif
#endif
