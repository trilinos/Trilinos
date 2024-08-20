// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
