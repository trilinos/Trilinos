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

/* This is similar to the Hash_i_dh class (woe, for a lack
   of templates); this this class is for hashing data
   consisting of single, non-negative integers.
*/


#ifndef HASH_I_DH
#define HASH_I_DH

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

/*
    class methods 
    note: all parameters are inputs; the only output 
          is the "int" returned by Hash_i_dhLookup.
*/
  extern void Hash_i_dhCreate (Hash_i_dh * h, int size);
  /* For proper operation, "size," which is the minimal
     size of the hash table, must be a power of 2.
     Or, pass "-1" to use the default.
   */


  extern void Hash_i_dhDestroy (Hash_i_dh h);
  extern void Hash_i_dhReset (Hash_i_dh h);

  extern void Hash_i_dhInsert (Hash_i_dh h, int key, int data);
  /* throws error if <data, data> is already inserted;
     grows hash table if out of space.
   */

  extern int Hash_i_dhLookup (Hash_i_dh h, int key);
  /* returns "data" associated with "key,"
     or -1 if "key" is not found.
   */

#ifdef __cplusplus
}
#endif
#endif
