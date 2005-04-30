/*
// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
//                 Copyright (2004) Sandia Corporation
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
// @HEADER
*/

#include "RTOp_parallel_helpers.h"

#define MY_MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MY_MAX(a,b) ( (a) > (b) ? (a) : (b) )

void RTOp_parallel_calc_overlap(
  RTOp_index_type global_dim_in, RTOp_index_type local_sub_dim_in, RTOp_index_type local_offset_in
  ,const RTOp_index_type first_ele_in, const RTOp_index_type sub_dim_in, const RTOp_index_type global_offset_in
  ,RTOp_index_type* overlap_first_local_ele, RTOp_index_type* overalap_local_sub_dim
  ,RTOp_index_type* overlap_global_offset
  )
{
  RTOp_index_type              global_sub_dim = 0;
#ifdef RTOp_DEBUG
  assert( overlap_first_local_ele );
  assert( overalap_local_sub_dim );
  assert( overlap_global_offset );
  /* ToDo: Check the rest of the preconditions! */
#endif
  /* Dimension of global sub-vector */
  global_sub_dim = sub_dim_in ? sub_dim_in : global_dim_in - (first_ele_in-1);
  /* */
  /* We need to determine if the local elements stored in this process overlap */
  /* with the global sub-vector that the client has requested. */
  /* */
  if( !( local_offset_in + local_sub_dim_in < first_ele_in
       || (first_ele_in-1) + global_sub_dim < local_offset_in + 1 ) )
  {
    /* */
    /* Determine how much of the local sub-vector stored in this process gets operated on. */
    /* If (first_ele_in-1) <= local_offset_in, then we start at the first element */
    /* in this process.  Otherwise, we need to to increment by first_ele_in - local_offset_in */
    /* */
    *overlap_first_local_ele = (first_ele_in-1) <= local_offset_in ? 1 : first_ele_in - local_offset_in;
    /* */
    /* Deterime the number of elements in the local sub-vector that overlap with the */
    /* requested logical sub-vector. */
    /* */
    *overalap_local_sub_dim  = (
      MY_MIN((first_ele_in-1)+global_sub_dim,local_offset_in+local_sub_dim_in) /* last overlap element in process */
      -
      MY_MAX(first_ele_in,local_offset_in+1)                                   /* first overlap element in process */
      + 1
      );
    /* */
    /* Finally, figure out where this local sub-vectors fit into the logical vector that the */
    /* client has specified with global_offset_in and first_ele_in.  Note that the element */
    /* this->(first_ele) acts as the the first element in the logical vector defined by the client */
    /* if gloabal_offset_in == 0.  Therefore, we need to subtract (first_ele_in - 1) from */
    /* local_offset_in to get the true offset into the logicl vector defined by the client.  Then */
    /* we can adjust it by adding global_offset_in to place it into the clients actual logical */
    /* vector.. */
    /* */
    *overlap_global_offset   = (
      ( first_ele_in - 1 > local_offset_in
        ? 0
        : local_offset_in - (first_ele_in - 1)
        )               /* First element in 'v' in logical sub-vector 'g' */
      + global_offset_in  /* Adding adjustment into logical sub-vector 'p' */
      );
  }
  else {
    *overlap_first_local_ele = 0; /* No overlap */
  }
}
