// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_RTOP_SUB_RANGE_DECORATOR_DEF_HPP
#define RTOPPACK_RTOP_SUB_RANGE_DECORATOR_DEF_HPP


#include "RTOpPack_RTOpSubRangeDecorator_decl.hpp"


namespace RTOpPack {


// Constructors, accessors


template<class Scalar>
RTOpSubRangeDecorator<Scalar>::RTOpSubRangeDecorator()
  : first_ele_offset_(0), sub_dim_(-1)
{}


template<class Scalar>
RTOpSubRangeDecorator<Scalar>::RTOpSubRangeDecorator(
  const RCP<RTOpT<Scalar> > &op,
  const Ordinal first_ele_offset,
  const Ordinal sub_dim
  )
  : first_ele_offset_(0), sub_dim_(-1)
{
  nonconstInitialize(op, first_ele_offset, sub_dim);
}


template<class Scalar>
RTOpSubRangeDecorator<Scalar>::RTOpSubRangeDecorator(
  const RCP<const RTOpT<Scalar> > &op,
  const Ordinal first_ele_offset,
  const Ordinal sub_dim
  )
  : first_ele_offset_(0), sub_dim_(-1)
{
  initialize(op, first_ele_offset, sub_dim);
}


template<class Scalar>
void
RTOpSubRangeDecorator<Scalar>::nonconstInitialize(
  const RCP<RTOpT<Scalar> > &op,
  const Ordinal first_ele_offset,
  const Ordinal sub_dim
  )
{
  op_.initialize(op);
  first_ele_offset_ = first_ele_offset;
  sub_dim_ = sub_dim;
}


template<class Scalar>
void
RTOpSubRangeDecorator<Scalar>::initialize(
  const RCP<const RTOpT<Scalar> > &op,
  const Ordinal first_ele_offset,
  const Ordinal sub_dim
  )
{
  op_.initialize(op);
  first_ele_offset_ = first_ele_offset;
  sub_dim_ = sub_dim;
}


template<class Scalar>
RCP<RTOpT<Scalar> >
RTOpSubRangeDecorator<Scalar>::getNonconstOp()
{
  return op_.getNonconstObj();
}


template<class Scalar>
RCP<const RTOpT<Scalar> >
RTOpSubRangeDecorator<Scalar>::getOp() const
{
  return op_.getConstObj();
}


// Overridden from RTOpT


template<class Scalar>
void RTOpSubRangeDecorator<Scalar>::get_reduct_type_num_entries_impl(
  const Ptr<int> &num_values,
  const Ptr<int> &num_indexes,
  const Ptr<int> &num_chars
  ) const
{
  op_->get_reduct_type_num_entries(num_values, num_indexes, num_chars);
}


template<class Scalar>
Teuchos::RCP<ReductTarget>
RTOpSubRangeDecorator<Scalar>::reduct_obj_create_impl() const
{
  return op_->reduct_obj_create();
}


template<class Scalar>
void RTOpSubRangeDecorator<Scalar>::reduce_reduct_objs_impl(
  const ReductTarget &in_reduct_obj, const Ptr<ReductTarget> &inout_reduct_obj
  ) const
{
  op_->reduce_reduct_objs(in_reduct_obj, inout_reduct_obj);
}


template<class Scalar>
void RTOpSubRangeDecorator<Scalar>::reduct_obj_reinit_impl(
  const Ptr<ReductTarget> &reduct_obj ) const
{
  op_->reduct_obj_reinit(reduct_obj);
}


template<class Scalar>
void RTOpSubRangeDecorator<Scalar>::extract_reduct_obj_state_impl(
  const ReductTarget &reduct_obj,
  const ArrayView<primitive_value_type> &value_data,
  const ArrayView<index_type> &index_data,
  const ArrayView<char_type> &char_data
  ) const
{
  op_->extract_reduct_obj_state(reduct_obj, value_data, index_data, char_data);
}


template<class Scalar>
void RTOpSubRangeDecorator<Scalar>::load_reduct_obj_state_impl(
  const ArrayView<const primitive_value_type> &value_data,
  const ArrayView<const index_type> &index_data,
  const ArrayView<const char_type> &char_data,
  const Ptr<ReductTarget> &reduct_obj
  ) const
{
  op_->load_reduct_obj_state(value_data, index_data, char_data, reduct_obj);
}


template<class Scalar>
std::string RTOpSubRangeDecorator<Scalar>::op_name_impl() const
{
  return (std::string("RTOpSubRangeDecorator{")+op_->op_name()+"}");
}


template<class Scalar>
bool RTOpSubRangeDecorator<Scalar>::coord_invariant_impl() const
{
  return op_->coord_invariant();
}


template<class Scalar>
void RTOpSubRangeDecorator<Scalar>::apply_op_impl(
  const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
  const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
  const Ptr<ReductTarget> &reduct_obj
  ) const
{
  
  // Check for full overlap
  if (first_ele_offset_ == 0 && sub_dim_ < 0) {
    // Entire range, just fall through
    op_->apply_op(sub_vecs, targ_sub_vecs, reduct_obj);
    return;
  }

  const Ordinal globalOffset =
    (sub_vecs.size() ? sub_vecs[0].globalOffset(): targ_sub_vecs[0].globalOffset());
  const Ordinal subDim = 
    (sub_vecs.size() ? sub_vecs[0].subDim(): targ_sub_vecs[0].subDim());

  // Check for no overlap
  if (globalOffset >= first_ele_offset_ + sub_dim_) {
    // No overlap
    return;
  }
  if (globalOffset + subDim <= first_ele_offset_) {
    // NO overlap
    return;
  }

  const Ordinal localOffset = 
    (first_ele_offset_ > globalOffset
      ? first_ele_offset_ - globalOffset
      : 0);

  const Ordinal localSubDim =
    std::min(globalOffset + subDim, first_ele_offset_ + sub_dim_)
    - (globalOffset + localOffset);
  
  Array<ConstSubVectorView<Scalar> > sub_sub_vecs(sub_vecs.size());
  for (int k = 0; k < sub_vecs.size(); ++k) {
    const Ordinal stride = sub_vecs[k].stride(); 
    sub_sub_vecs[k].initialize(
      globalOffset+ localOffset,
      localSubDim,
      sub_vecs[k].values().persistingView(localOffset*stride, localSubDim*stride),
      stride
      );
  }

  Array<SubVectorView<Scalar> > targ_sub_sub_vecs(targ_sub_vecs.size());
  for (int k = 0; k < targ_sub_vecs.size(); ++k) {
    const Ordinal stride = targ_sub_vecs[k].stride(); 
    targ_sub_sub_vecs[k].initialize(
      globalOffset+ localOffset,
      localSubDim,
      targ_sub_vecs[k].values().persistingView(localOffset*stride, localSubDim*stride),
      stride
      );
  }

  op_->apply_op(sub_sub_vecs(), targ_sub_sub_vecs(), reduct_obj);

}


} // namespace RTOpPack


#endif // RTOPPACK_RTOP_SUB_RANGE_DECORATOR_DEF_HPP
