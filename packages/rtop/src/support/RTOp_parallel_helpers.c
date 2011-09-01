/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "RTOp_parallel_helpers.h"

#define MY_MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MY_MAX(a,b) ( (a) > (b) ? (a) : (b) )

void RTOp_parallel_calc_overlap(
  Teuchos_Index global_dim_in, Teuchos_Index local_sub_dim_in, Teuchos_Index local_off_in
  ,const Teuchos_Index first_ele_off_in, const Teuchos_Index sub_dim_in, const Teuchos_Index global_off_in
  ,Teuchos_Index* overlap_first_local_ele_off, Teuchos_Index* overalap_local_sub_dim
  ,Teuchos_Index* overlap_global_off
  )
{
  Teuchos_Index  global_sub_dim = 0;
#ifdef RTOp_DEBUG
  assert( overlap_first_local_ele_off );
  assert( overalap_local_sub_dim );
  assert( overlap_global_off );
  /* ToDo: Check the rest of the preconditions! */
#endif
  /* Dimension of global sub-vector */
  global_sub_dim = sub_dim_in >= 0 ? sub_dim_in : global_dim_in - first_ele_off_in;
  /*
   * We need to determine if the local elements stored in this process overlap
   * with the global sub-vector that the client has requested.
   */
  if( !(
        local_off_in + local_sub_dim_in < first_ele_off_in + 1
        ||
        first_ele_off_in + global_sub_dim < local_off_in + 1
        )
    )
  {
    /* 
     * Determine how much of the local sub-vector stored in this process gets operated on. 
     * If (first_ele_off_in-1) <= local_off_in, then we start at the first element 
     * in this process.  Otherwise, we need to to increment by first_ele_off_in - local_off_in 
     */
    *overlap_first_local_ele_off = first_ele_off_in <= local_off_in ? 0 : first_ele_off_in - local_off_in;
    /*
     * Deterime the number of elements in the local sub-vector that overlap with the
     * requested logical sub-vector.
     */
    *overalap_local_sub_dim  = (
      MY_MIN(first_ele_off_in+global_sub_dim,local_off_in+local_sub_dim_in) /* last overlap element plus 1 in process */
      -
      MY_MAX(first_ele_off_in,local_off_in)                                 /* first overlap element in process */
      );
    /*
     * Finally, figure out where this local sub-vectors fit into the logical
     * vector that the client has specified with global_off_in and
     * first_ele_off_in.  Note that the element this->(first_ele_off) acts as
     * the the first element in the logical vector defined by the client if
     * gloabal_offset_in == 0.  Therefore, we need to subtract
     * first_ele_off_in from local_off_in to get the true offset into the
     * logicl vector defined by the client.  Then we can adjust it by adding
     * global_off_in to place it into the clients actual logical vector..
     */
    *overlap_global_off   = (
      ( first_ele_off_in > local_off_in
        ? 0
        : local_off_in - first_ele_off_in
        )                 /* First element in 'v' in logical sub-vector 'g' */
      + global_off_in     /* Adding adjustment into logical sub-vector 'p' */
      );
  }
  else {
    *overlap_first_local_ele_off = -1; /* No overlap */
  }
}
