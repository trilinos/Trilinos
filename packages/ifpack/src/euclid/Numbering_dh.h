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

#ifndef NUMBERING_DH_H
#define NUMBERING_DH_H


/* code and algorithms in this class adopted from Edmond Chow's
   ParaSails
*/


#include "euclid_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _numbering_dh {
  int   size;    /* max number of indices that can be stored;
                    (length of idx_ext[]) 
                  */
  int   first;   /* global number of 1st local index (row) */
  int   m;       /* number of local indices (number of local rows in mat) */
  int   *idx_ext;   /* sorted list of external indices */
  int   *idx_extLo; /* sorted list of external indices that are < first */
  int   *idx_extHi; /* sorted list of external indices that are >= first+m */
  int   num_ext; /* number of external (non-local) indices = num_extLo+num_extHi */
  int   num_extLo; /* number of external indices < first */
  int   num_extHi; /* number of external indices >= first+num_loc */
  Hash_i_dh global_to_local;

  bool debug;
};

extern void Numbering_dhCreate(Numbering_dh *numb);
extern void Numbering_dhDestroy(Numbering_dh numb);

  /* must be called before calling Numbering_dhGlobalToLocal() or
     Numbering_dhLocalToGlobal().
   */
extern void Numbering_dhSetup(Numbering_dh numb, Mat_dh mat);


  /* input: global_in[len], which contains global row numbers.
     output: local_out[len], containing corresponding local numbers.
     note: global_in[] and local_out[] may be identical.
   */
extern void Numbering_dhGlobalToLocal(Numbering_dh numb, int len, 
                                      int *global_in, int *local_out);

#ifdef __cplusplus
}
#endif
#endif
