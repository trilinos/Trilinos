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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_MPI_APPLY_OP_DECL_HPP
#define RTOPPACK_MPI_APPLY_OP_DECL_HPP

#include "RTOpPack_RTOpT.hpp"
#include "RTOp_MPI_config.h"

#define RTOPPACK_MPI_APPLY_OP_DUMP

namespace RTOpPack {

#ifdef RTOPPACK_MPI_APPLY_OP_DUMP
extern bool show_mpi_apply_op_dump;
#endif // RTOPPACK_MPI_APPLY_OP_DUMP

/** \brief Initialize MPI compatible type signature arrays for
 * reduction/transformation operator object instance data and
 * reduction target object data.
 *
 * @param num_values    [in] Number of primitive_value members
 * @param num_indexes   [in] Number of index members
 * @param num_chars     [in] Number of character members
 * @param num_entries   [out] Number of entries in output arrays set
 * @param block_lengths [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 * @param displacements [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 * @param datatypes     [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 *
 * See the MPI function <tt>MPI_Type_struct(...)</tt> for a description of these arrays.
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class primitive_value_type>
void MPI_type_signature(
  const int num_values
  ,const int num_indexes
  ,const int num_chars
  ,int* num_entries
  ,int block_lengths[]
  ,MPI_Aint displacements[]
  ,MPI_Datatype datatypes[]
  );

/** \brief Return the size in bytes of an external representation of <tt>reduct_obj</tt>.
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class primitive_value_type>
int reduct_obj_ext_size(
  int   num_values
  ,int  num_indexes
  ,int  num_chars
  )
{
  return (3 + num_values) * sizeof(primitive_value_type)
    + num_indexes         * sizeof(index_type)
    + num_chars           * sizeof(char_type);
}

/** \brief .
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void extract_reduct_obj_ext_state(
  const RTOpT<Scalar>    &op
  ,const ReductTarget    &reduct_obj
  ,int                   num_values
  ,int                   num_indexes
  ,int                   num_chars
  ,void                  *reduct_obj_ext
  );

/** \brief .
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void load_reduct_obj_ext_state(
  const RTOpT<Scalar>    &op
  ,const void            *reduct_obj_ext
  ,ReductTarget          *reduct_obj
  );

/** \brief Apply an RTOp in SMPD mode with MPI to a set of vectors with
 * contiguous storage per processor.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void MPI_apply_op(
  MPI_Comm                            comm
  ,const RTOpT<Scalar>                &op
  ,const int                          root_rank
  ,const int                          num_vecs
  ,const ConstSubVectorView<Scalar>   sub_vecs[]
  ,const int                          num_targ_vecs
  ,const SubVectorView<Scalar>        targ_sub_vecs[]
  ,ReductTarget                       *reduct_obj
  );

/** \brief Apply an RTOp in SMPD mode with MPI to a set of columns to a set of
 * multi-vectors with contiguous storage per processor.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void MPI_apply_op(
  MPI_Comm                                 comm
  ,const RTOpT<Scalar>                     &op
  ,const int                               root_rank
  ,const int                               num_cols
  ,const int                               num_multi_vecs
  ,const ConstSubMultiVectorView<Scalar>   sub_multi_vecs[]
  ,const int                               num_targ_multi_vecs
  ,const SubMultiVectorView<Scalar>        targ_sub_multi_vecs[]
  ,ReductTarget*const                      reduct_objs[]
  );

/** \brief Perform global reduction of reduction target objects.
 *
 * \param  comm   [in] MPI communicator
 * \param  op     [in] Reduction operator object associated with reduction objects.
 *                This operator defines how the global reductions are performed.
 * \param  root_rank
 *                [in] The rank of the root processor
 * \param  num_cols
 *                [in] The number of vector sets (i.e. multi-vector columns) that these
 *                reduction objects where collected over.
 * \param  i_reduct_objs
 *                [in] Array (length <tt>num_cols</tt>) of local intermediate reduction objects.
 * \param  reduct_objs
 *                [in/out] Array (length <tt>num_cols</tt>) of global reduction objects
 *                that each <tt>i_reduct_objs[i]<tt> from each processor will be reduced into.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void MPI_all_reduce(
  MPI_Comm                            comm
  ,const RTOpT<Scalar>                &op
  ,const int                          root_rank
  ,const int                          num_cols
  ,const ReductTarget*const           i_reduct_objs[]
  ,ReductTarget*const                 reduct_objs[]
  );

/** \brief Apply an RTOp in SMPD mode with MPI to a set of columns to a set of
 * multi-vectors with contiguous storage per processor.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void  MPI_apply_op(
  MPI_Comm                                  comm
  ,const RTOpT<Scalar>                      &op
  ,const int                                root_rank
  ,const int                                num_cols
  ,const int                                num_vecs
  ,const ConstSubVectorView<Scalar>         sub_vecs[]
  ,const int                                num_targ_vecs
  ,const SubVectorView<Scalar>              sub_targ_vecs[]
  ,ReductTarget*const                       reduct_objs[]
  );

} // end namespace RTOpPack

#endif // RTOPPACK_MPI_APPLY_OP_DECL_HPP
