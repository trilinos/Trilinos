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

#ifndef NUMBERING_DH_H
#define NUMBERING_DH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif


/* code and algorithms in this class adopted from Edmond Chow's
   ParaSails
*/


#include "euclid_common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  struct _numbering_dh
  {
    int size;			/* max number of indices that can be stored;
				   (length of idx_ext[]) 
				 */
    int first;			/* global number of 1st local index (row) */
    int m;			/* number of local indices (number of local rows in mat) */
    int *idx_ext;		/* sorted list of external indices */
    int *idx_extLo;		/* sorted list of external indices that are < first */
    int *idx_extHi;		/* sorted list of external indices that are >= first+m */
    int num_ext;		/* number of external (non-local) indices = num_extLo+num_extHi */
    int num_extLo;		/* number of external indices < first */
    int num_extHi;		/* number of external indices >= first+num_loc */
    Hash_i_dh global_to_local;

    bool debug;
  };

  extern void Numbering_dhCreate (Numbering_dh * numb);
  extern void Numbering_dhDestroy (Numbering_dh numb);

  /* must be called before calling Numbering_dhGlobalToLocal() or
     Numbering_dhLocalToGlobal().
   */
  extern void Numbering_dhSetup (Numbering_dh numb, Mat_dh mat);


  /* input: global_in[len], which contains global row numbers.
     output: local_out[len], containing corresponding local numbers.
     note: global_in[] and local_out[] may be identical.
   */
  extern void Numbering_dhGlobalToLocal (Numbering_dh numb, int len,
					 int *global_in, int *local_out);

#ifdef __cplusplus
}
#endif
#endif
