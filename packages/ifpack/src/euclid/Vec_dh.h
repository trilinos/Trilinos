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

#ifndef VEC_DH_H
#define VEC_DH_H

#include "euclid_common.h"
#ifdef __cplusplus
extern "C" {
#endif

struct _vec_dh {
  int n;
  double *vals;
};

extern void Vec_dhCreate(Vec_dh *v);
extern void Vec_dhDestroy(Vec_dh v);
extern void Vec_dhInit(Vec_dh v, int size);
        /* allocates storage, but does not initialize values */

extern void Vec_dhDuplicate(Vec_dh v, Vec_dh *out);
        /* creates vec and allocates storage, but neither
         * initializes nor copies values 
         */

extern void Vec_dhCopy(Vec_dh x, Vec_dh y);
        /* copies values from x to y;
         * y must have proper storage allocated,
         * e.g, through previous call to Vec_dhDuplicate,
         * or Vec_dhCreate and Vec_dhInit.
         */

extern void Vec_dhSet(Vec_dh v, double value);
extern void Vec_dhSetRand(Vec_dh v);

extern void Vec_dhRead(Vec_dh *v, int ignore, char *filename);
extern void Vec_dhReadBIN(Vec_dh *v, char *filename);
extern void Vec_dhPrint(Vec_dh v, SubdomainGraph_dh sg, char *filename);
extern void Vec_dhPrintBIN(Vec_dh v, SubdomainGraph_dh sg, char *filename); 
#ifdef __cplusplus
}
#endif
#endif
