// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_APPLY_OP_HELPER_DECL_HPP
#define THYRA_APPLY_OP_HELPER_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "RTOpPack_RTOpT.hpp"


namespace Thyra {


/** \brief Validate the inputs to <tt>VectorBase::applyOp()</tt>.
 *
 * Throws an exception with a nice error message if one of the
 * preconditions are not met.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
void apply_op_validate_input(
  const std::string &func_name,
  const VectorSpaceBase<Scalar> &space,
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset
  );


/** \brief Validate the inputs to <tt>MultiVectorBase::applyOp()</tt>.
 *
 * Throws an exception with a nice error message if one of the
 * preconditions are not met.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
void apply_op_validate_input(
  const std::string &func_name,
  const VectorSpaceBase<Scalar> &domain,
  const VectorSpaceBase<Scalar> &range,
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal primary_global_offset
  );


} // end namespace Thyra


#endif // THYRA_APPLY_OP_HELPER_DECL_HPP
