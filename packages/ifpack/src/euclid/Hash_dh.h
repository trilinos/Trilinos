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

#ifndef HASH_D_DH
#define HASH_D_DH

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
extern "C" {
#endif

typedef struct _hash_node {
  int     iData;      /* integer */
  double  fData;      /* float */
  int     *iDataPtr;  /* pointer to integer */
  int     *iDataPtr2; /* pointer to integer */
  double  *fDataPtr;  /* pointer to float */
} HashData;


typedef struct _hash_node_private HashRecord;

/* data structure for the hash table; do not directly access */
struct _hash_dh {
  int         size;   /* total slots in table */
  int         count;  /* number of insertions in table */
  int         curMark;
  HashRecord  *data;
};


extern void Hash_dhCreate(Hash_dh *h, int size);
extern void Hash_dhDestroy(Hash_dh h);
extern void Hash_dhInsert(Hash_dh h, int key, HashData *data);
extern HashData *Hash_dhLookup(Hash_dh h, int key);
         /* returns NULL if record isn't found. */

extern void Hash_dhReset(Hash_dh h);
extern void Hash_dhPrint(Hash_dh h, FILE *fp);


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
