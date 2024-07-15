// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef RTOP_PARALLEL_HELPERS_H
#define RTOP_PARALLEL_HELPERS_H

#include "RTOp_ConfigDefs.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief This function helps to implement vector method
 * <tt>apply_op(...)</tt> for any type of parallel vector.
 *
 * \param global_dim [in] Dimension of the original parallel vector 'v' (see
 * above).
 *
 * \param local_sub_dim [in] Dimension of the local subvector 'u' (see above).
 *
 * \param local_off [in] Gives the offset of the first element in the local
 * sub-vector 'u' into the global vector 'v' (see above).
 *
 * \param first_ele_off [in] Determines the first element in 'v' which is used
 * to define the logical sub-vector 'g' (see above).
 *
 * \param sub_dim [in] Determines the length of the logical sub-vector 'g'
 * (see above).  If <tt>sub_dim < 0</tt> then <tt>sub_dim = global_dim -
 * first_ele_off</tt> is used in its place.
 *
 * \param global_off [in] Determines the offset of the logical subvector 'g'
 * into the logical global vector 'p' (see above).
 *
 * \param overlap_first_local_ele_off [out] If <tt>*overlap_first_local_ele <
 * 0</tt> on output, then this means that there is no overlap of 'u' with 'g'
 * (see above).  Otherwise, there is overlap and
 * <tt>*overlap_first_local_ele_off</tt> gives the first element in 'u' that
 * overlaps with 'g' which defines 'w' (see above).
 *
 * \param overlap_local_sub_dim [out] If <tt>*overlap_first_local_ele_off <
 * 0</tt> on output then this argument is not set and should be ignored.
 * Otherwise, <tt>*overlap_local_sub_dim</tt> gives number of elements in 'u'
 * that overlaps with 'g' that defines 'w' (see above).
 *
 * \param overlap_global_off [out] If <tt>*overlap_first_local_ele_off <
 * 0</tt> on output then this argument is not set and should be ignored.
 * Otherwise, <tt>*overlap_global_off</tt> gives the placement of 'w' into 'p'
 * (see above).
 *
 * Preconditions:<ul>
 * <li><tt>global_dim > 0</tt>
 * <li><tt>local_sub_dim > 0</tt>
 * <li><tt>0 <= local_off <= global_dim - local_sub_dim</tt>
 * <li><tt>1 <= first_ele_off <= global_dim</tt>
 * <li>[<tt>sub_dim > 0</tt>] <tt>0 < sub_dim <= global_dim - first_ele_off</tt>
 * <li><tt>0 <= global_off</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>0 <= overlap_first_local_ele_off <= local_sub_dim</tt>
 * <li>[<tt>overlap_first_local_ele_off == 0</tt>] There is no overlap of 'g' with 'u'
 * <li>[<tt>overlap_first_local_ele_off != 0</tt>] <tt>0 <= overlap_local_sub_dim
 *   <= local_sub_dim - overlap_first_local_ele_off</tt>
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
 * input arguments \c first_ele_off, \c sub_dim and \c global_off.  If for the time being
 * we ignore \c global_off, and then 'g' is defined in terms of 'v' as:
 \verbatim
 
 g(k) = v(first_ele_off+k), for k = 0...sub_dim-1
 \endverbatim
 * However, for greater flexibility, the client can specify that the logical vector 'g'
 * is really a sub-vector in a larger vector (a third vector, call it 'p') and can therefore
 * specify where 'g' exists in 'p' using \c global_off as:
 \verbatim
 
 p(k+global_off) = g(k) = v(first_ele_off+k), for k = 0...sub_dim-1
 \endverbatim
 * In order to apply a reduction/transformation operator over the sub-vector 'g' in 'p'
 * each process can only work with the elements of 'v' stored in the local process.  Specifically,
 * the local elements of 'v' stored in this process (the fourth vector, all it 'u') are:
 \verbatim
 
 u(k) = v(local_off+k), for k = 0...local_sub_dim-1
 \endverbatim
 * The tricky part of implementing this function is is determining how much of 'u' overlaps
 * with 'g' and then getting the offset into 'p' correct.  If the local elements 'u' overlaps
 * with 'g' then this defines the fifth sub-vector (call it 'w') that defines the overlap
 * and is specified by the return arguments as:
 \verbatim
 
 w(k) = p(overlap_global_off+k) = u(overlap_first_local_ele_off+k), for k = 0...overalap_local_sub_dim-1
 \endverbatim
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
void RTOp_parallel_calc_overlap(
  Teuchos_Ordinal global_dim, Teuchos_Ordinal local_sub_dim, Teuchos_Ordinal local_off
  ,const Teuchos_Ordinal first_ele_off, const Teuchos_Ordinal sub_dim, const Teuchos_Ordinal global_off
  ,Teuchos_Ordinal* overlap_first_local_ele_off, Teuchos_Ordinal* overalap_local_sub_dim
  ,Teuchos_Ordinal* overlap_global_off
  );

#ifdef __cplusplus
}
#endif

#endif /* RTOP_PARALLEL_HELPERS_H */
