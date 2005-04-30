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

#ifndef RTOP_PARALLEL_HELPERS_H
#define RTOP_PARALLEL_HELPERS_H

#include "RTOp_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief This function helps to implement vector method <tt>apply_op(...)</tt> for any type of parallel vector.
 *
 * @param  global_dim    [in] Dimension of the original parallel vector 'v' (see above). 
 * @param  local_sub_dim [in] Dimension of the local subvector 'u' (see above).
 * @param  local_offset  [in] Gives the offset of the first element in the local sub-vector 'u'
 *                       into the global vector 'v' (see above).
 * @param  first_ele     [in] Determines the first element in 'v' which is used to define
 *                       the logical sub-vector 'g' (see above).
 * @param  sub_dim       [in] Determines the length of the logical sub-vector 'g' (see above).
 *                       If <tt>sub_dim == 0</tt> then <tt>sub_dim = global_dim - (first_ele-1)</tt>
 *                       is used in its place.
 * @param  global_offset [in] Determines the offset of the logical subvector 'g' into the logical
 *                       global vector 'p' (see above).
 * @param  overlap_first_local_ele
 *                       [out] If <tt>*overlap_first_local_ele == 0</tt> on output, then this means
 *                       that there is no overlap of 'u' with 'g' (see above).  Otherwise, there
 *                       is overlap and <tt>*overlap_first_local_ele</tt> gives the first element
 *                       in 'u' that overlaps with 'g' which defines 'w' (see above).
 * @param  overlap_local_sub_dim
 *                       [out]  If <tt>*overlap_first_local_ele == 0</tt> on output then this
 *                       argument is not set and should be ignored.  Otherwise, 
 *                       <tt>*overlap_local_sub_dim</tt> gives number of elements in 'u' that
 *                       overlaps with 'g' that defines 'w' (see above).
 * @param  overlap_global_offset
 *                       [out]  If <tt>*overlap_first_local_ele == 0</tt> on output then this
 *                       argument is not set and should be ignored.  Otherwise, 
 *                       <tt>*overlap_global_offset</tt> gives the placement of 'w' into 'p'
 *                       (see above).
 *
 * Preconditions:<ul>
 * <li><tt>global_dim > 0</tt>
 * <li><tt>local_sub_dim > 0</tt>
 * <li><tt>0 <= local_offset <= global_dim - local_sub_dim</tt>
 * <li><tt>1 <= first_ele <= global_dim</tt>
 * <li>[<tt>sub_dim != 0</tt>] <tt>0 < sub_dim <= global_dim - (first_ele-1)</tt>
 * <li><tt>0 <= global_offset</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>0 <= overlap_first_local_ele <= local_sub_dim</tt>
 * <li>[<tt>overlap_first_local_ele == 0</tt>] There is no overlap of 'g' with 'u'
 * <li>[<tt>overlap_first_local_ele != 0</tt>] <tt>0 <= overlap_local_sub_dim <= local_sub_dim - (overlap_first_local_ele-1)</tt>
 * </ul>
 *
 * To understand what this function computes first consider the what
 * an <tt>apply_op(...)</tt> method might look like from a vector
 * class. This method would take a list of non-mutable vectors:

 \verbatim

 v_1, v_2, ..., v_p
 \endverbatim
 * and a list of mutable vectors
 \verbatim

 z_1, z_2, ..., z_q
 \endverbatim
 * and then apply the reduction/transformation operator \c op over some subset
 * of the elements in these vectors according to their placement as a set of
 * sub-vectors in some other logical vector.
 * 
 * Let's consider how things are treated for a single vector argument \c v_i or \c z_i
 * which we will call 'v'.  This global vector 'v' is the first vector that we identity.
 * One must understand that there are five
 * distict vectors (or sub-vectors) being refered to here.  The first vector (call it 'v')
 * and is one of the parallel vectors that <tt>apply_op()</tt> is called on that we have
 * already discussed.  The second vector (call it 'g') is the
 * logical sub-vector that the client wants to
 * represent using the elements in 'v'.  This logical sub-vector is specified by the
 * input arguments \c first_ele, \c sub_dim and \c global_offset.  If for the time being
 * we ignore \c global_offset, and then 'g' is defined in terms of 'v' as:
 \verbatim
 
 g(k) = v((first_ele-1)+k), for k = 1...sub_dim
 \endverbatim
 * However, for greater flexibility, the client can specify that the logical vector 'g'
 * is really a sub-vector in a larger vector (a third vector, call it 'p') and can therefore
 * specify where 'g' exists in 'p' using \c global_offset as:
 \verbatim
 
 p(k+global_offset) = g(k) = v((first_ele-1)+k), for k = 1...sub_dim
 \endverbatim
 * In order to apply a reduction/transformation operator over the sub-vector 'g' in 'p'
 * each process can only work with the elements of 'v' stored in the local process.  Specifically,
 * the local elements of 'v' stored in this process (the fourth vector, all it 'u') are:
 \verbatim
 
 u(k) = v(local_offset+k), for k = 1...local_sub_dim
 \endverbatim
 * The tricky part of implementing this function is is determining how much of 'u' overlaps
 * with 'g' and then getting the offset into 'p' correct.  If the local elements 'u' overlaps
 * with 'g' then this defines the fifth sub-vector (call it 'w') that defines the overlap
 * and is specified by the return arguments as:
 \verbatim
 
 w(k) = p(overlap_global_offset+k) = u((overlap_first_local_ele-1)+k), for k = 1...overalap_local_sub_dim
 \endverbatim
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
void RTOp_parallel_calc_overlap(
  RTOp_index_type global_dim, RTOp_index_type local_sub_dim, RTOp_index_type local_offset
  ,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
  ,RTOp_index_type* overlap_first_local_ele, RTOp_index_type* overalap_local_sub_dim
  ,RTOp_index_type* overlap_global_offset
  );

#ifdef __cplusplus
}
#endif

#endif /* RTOP_PARALLEL_HELPERS_H */
