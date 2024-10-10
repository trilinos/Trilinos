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

#ifndef HASH_D_DH
#define HASH_D_DH

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

/* todo: rehash should be implemented in Hash_dhInsert();
         as of now, an error is set if the table overflows.
*/

#include "euclid_common.h"

/* This should be done with templates, if this were in C++;
   for now, a record contains every type of entry we might
   need; this is a waste of memory, when one is only intersted
   in hashing <key, int> pairs!
*/
#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct _hash_node
  {
    int iData;			/* integer */
    double fData;		/* float */
    int *iDataPtr;		/* pointer to integer */
    int *iDataPtr2;		/* pointer to integer */
    double *fDataPtr;		/* pointer to float */
  } HashData;


  typedef struct _hash_node_private HashRecord;

/* data structure for the hash table; do not directly access */
  struct _hash_dh
  {
    int size;			/* total slots in table */
    int count;			/* number of insertions in table */
    int curMark;
    HashRecord *data;
  };


  extern void Hash_dhCreate (Hash_dh * h, int size);
  extern void Hash_dhDestroy (Hash_dh h);
  extern void Hash_dhInsert (Hash_dh h, int key, HashData * data);
  extern HashData *Hash_dhLookup (Hash_dh h, int key);
  /* returns NULL if record isn't found. */

  extern void Hash_dhReset (Hash_dh h);
  extern void Hash_dhPrint (Hash_dh h, FILE * fp);


#define HASH_1(k,size,idxOut)    \
         {  *idxOut = k % size;  }

#define HASH_2(k,size,idxOut)      \
          {  \
            int r = k % (size-13); \
            r = (r % 2) ? r : r+1; \
            *idxOut = r;           \
          }

#ifdef __cplusplus
}
#endif
#endif
