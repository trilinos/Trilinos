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

/* This is similar to the Hash_i_dh class (woe, for a lack
   of templates); this this class is for hashing data
   consisting of single, non-negative integers.
*/


#ifndef HASH_I_DH
#define HASH_I_DH

#include "euclid_common.h"

                                 
#ifdef __cplusplus
extern "C" {
#endif

/*
    class methods 
    note: all parameters are inputs; the only output 
          is the "int" returned by Hash_i_dhLookup.
*/
extern void Hash_i_dhCreate(Hash_i_dh *h, int size);
  /* For proper operation, "size," which is the minimal
     size of the hash table, must be a power of 2.
     Or, pass "-1" to use the default.
   */


extern void Hash_i_dhDestroy(Hash_i_dh h);
extern void Hash_i_dhReset(Hash_i_dh h);

extern void Hash_i_dhInsert(Hash_i_dh h, int key, int data);
  /* throws error if <data, data> is already inserted;
     grows hash table if out of space.
   */

extern int  Hash_i_dhLookup(Hash_i_dh h, int key);
    /* returns "data" associated with "key,"
       or -1 if "key" is not found.
     */

#ifdef __cplusplus
}
#endif
#endif
