/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef VEC_DH_H
#define VEC_DH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "euclid_common.h"
#ifdef __cplusplus
extern "C"
{
#endif

  struct _vec_dh
  {
    int n;
    double *vals;
  };

  extern void Vec_dhCreate (Vec_dh * v);
  extern void Vec_dhDestroy (Vec_dh v);
  extern void Vec_dhInit (Vec_dh v, int size);
  /* allocates storage, but does not initialize values */

  extern void Vec_dhDuplicate (Vec_dh v, Vec_dh * out);
  /* creates vec and allocates storage, but neither
   * initializes nor copies values 
   */

  extern void Vec_dhCopy (Vec_dh x, Vec_dh y);
  /* copies values from x to y;
   * y must have proper storage allocated,
   * e.g, through previous call to Vec_dhDuplicate,
   * or Vec_dhCreate and Vec_dhInit.
   */

  extern void Vec_dhSet (Vec_dh v, double value);
  extern void Vec_dhSetRand (Vec_dh v);

  extern void Vec_dhRead (Vec_dh * v, int ignore, char *filename);
  extern void Vec_dhReadBIN (Vec_dh * v, char *filename);
  extern void Vec_dhPrint (Vec_dh v, SubdomainGraph_dh sg, char *filename);
  extern void Vec_dhPrintBIN (Vec_dh v, SubdomainGraph_dh sg, char *filename);
#ifdef __cplusplus
}
#endif
#endif
