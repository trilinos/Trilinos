// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_GET_SUB_VECTOR_DEF_HPP
#define RTOPPACK_ROP_GET_SUB_VECTOR_DEF_HPP


#include "RTOpPack_ROpGetSubVector_decl.hpp"


namespace RTOpPack {


template<class Scalar>
ROpGetSubVector<Scalar>::ROpGetSubVector( const index_type l,
  const index_type u
  )
  :RTOpT<Scalar>("ROpGetSubVector"), l_(l), u_(u)
{}


template<class Scalar>
void ROpGetSubVector<Scalar>::set_range( const index_type l,
  const index_type u
  )
{
  l_ = l;
  u_ = u;
}


template<class Scalar>
const ConstSubVectorView<Scalar>
ROpGetSubVector<Scalar>::operator()( const ReductTarget& reduct_obj ) const
{
  using Teuchos::dyn_cast;
  return dyn_cast<const DefaultReductTarget<SubVectorView< Scalar> > >(reduct_obj).get();
}


// Overridden from RTOpT


template<class Scalar>
void ROpGetSubVector<Scalar>::get_reduct_type_num_entries_impl(
  const Ptr<int> &num_values,
  const Ptr<int> &num_indexes,
  const Ptr<int> &num_chars
  ) const
{
  typedef PrimitiveTypeTraits<Scalar,Scalar> PTT;
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  *num_values = (u_-l_+1)*num_prim_objs_per_scalar;
  *num_indexes = 0;
  *num_chars = 0;
}


template<class Scalar>
Teuchos::RCP<ReductTarget>
ROpGetSubVector<Scalar>::reduct_obj_create_impl() const
{
  const index_type subDim = u_ - l_ + 1;
  const ArrayRCP<Scalar> values = Teuchos::arcp<Scalar>(subDim);
  std::fill(values.begin(), values.end(), ScalarTraits<Scalar>::zero());
  return defaultReductTarget(
    SubVectorView<Scalar>( l_, subDim, values, 1 )
    );
}


template<class Scalar>
void ROpGetSubVector<Scalar>::reduce_reduct_objs_impl(
  const ReductTarget &in_reduct_obj, const Ptr<ReductTarget> &inout_reduct_obj
  ) const
{

  using Teuchos::dyn_cast;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;

  DRTSVV &drtsvv_inout_reduct_obj = dyn_cast<DRTSVV>(*inout_reduct_obj);

  const ConstSubVectorView<Scalar> sub_vec_in =
    dyn_cast<const DRTSVV>(in_reduct_obj).get();
  SubVectorView<Scalar> sub_vec_inout = drtsvv_inout_reduct_obj.get();

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(
    sub_vec_in.subDim()!=sub_vec_inout.subDim()
    || sub_vec_in.globalOffset()!=sub_vec_inout.globalOffset()
    || is_null(sub_vec_in.values())
    || is_null(sub_vec_inout.values())
    || sub_vec_in.stride()!=1
    || sub_vec_inout.stride()!=1
    );
#endif // TEUCHOS_DEBUG

  typedef typename ArrayRCP<const Scalar>::const_iterator const_iter_t;
  typedef typename ArrayRCP<Scalar>::iterator iter_t;

  const_iter_t in_iter = sub_vec_in.values().begin();
  iter_t inout_iter = sub_vec_inout.values().begin();

  for( int k = 0; k < sub_vec_in.subDim(); ++k ) {
    *inout_iter++ += *in_iter++;
  }

  drtsvv_inout_reduct_obj.set(sub_vec_inout);

}


template<class Scalar>
void ROpGetSubVector<Scalar>::reduct_obj_reinit_impl(
  const Ptr<ReductTarget> &reduct_obj ) const
{
  using Teuchos::dyn_cast;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;
  DRTSVV &drtsvv_inout_reduct_obj = dyn_cast<DRTSVV>(*reduct_obj);
  SubVectorView<Scalar> sub_vec = drtsvv_inout_reduct_obj.get();
  std::fill( sub_vec.values().begin(), sub_vec.values().end(),
    ScalarTraits<Scalar>::zero() );
}


template<class Scalar>
void ROpGetSubVector<Scalar>::extract_reduct_obj_state_impl(
  const ReductTarget &reduct_obj,
  const ArrayView<primitive_value_type> &value_data,
  const ArrayView<index_type> &/* index_data */,
  const ArrayView<char_type> &/* char_data */
  ) const
{
  using Teuchos::null;
  using Teuchos::dyn_cast;
  typedef PrimitiveTypeTraits<Scalar,Scalar> PTT;
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  const ConstSubVectorView<Scalar> sub_vec =
    dyn_cast<const DefaultReductTarget<SubVectorView<Scalar> > >(reduct_obj).get();
  int value_data_off = 0;
  for(
    int k = 0;
    k < sub_vec.subDim();
    ++k, value_data_off += num_prim_objs_per_scalar
    )
  {
    PTT::extractPrimitiveObjs( sub_vec[k],
      value_data(value_data_off, num_prim_objs_per_scalar),
      null, null );
  }
}


template<class Scalar>
void ROpGetSubVector<Scalar>::load_reduct_obj_state_impl(
  const ArrayView<const primitive_value_type> &value_data,
  const ArrayView<const index_type> &/* index_data */,
  const ArrayView<const char_type> &/* char_data */,
  const Ptr<ReductTarget> &reduct_obj
  ) const
{
  using Teuchos::null;
  using Teuchos::outArg;
  using Teuchos::dyn_cast;
  using Teuchos::arcp_const_cast;
  typedef PrimitiveTypeTraits<Scalar,Scalar> PTT;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  DRTSVV &drtsvv_reduct_obj = dyn_cast<DRTSVV>(*reduct_obj);
  const ConstSubVectorView<Scalar> const_sub_vec = drtsvv_reduct_obj.get();
  const ArrayRCP<Scalar> sv_values =
    arcp_const_cast<Scalar>(const_sub_vec.values());
  int value_data_off = 0;
  for(
    int k = 0;
    k < const_sub_vec.subDim();
    ++k, value_data_off += num_prim_objs_per_scalar
    )
  {
    PTT::loadPrimitiveObjs(
      value_data(value_data_off, num_prim_objs_per_scalar), null, null,
      outArg(sv_values[k]) );
  }
}


template<class Scalar>
bool ROpGetSubVector<Scalar>::coord_invariant_impl() const
{
  return false;
}


template<class Scalar>
void ROpGetSubVector<Scalar>::apply_op_impl(
  const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
  const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
  const Ptr<ReductTarget> &reduct_obj
  ) const
{

  using Teuchos::dyn_cast;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;

  validate_apply_op( *this, 1, 0, true,
    sub_vecs, targ_sub_vecs, reduct_obj.getConst() );

  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
  const index_type subDim  = sub_vecs[0].subDim();
  const index_type globalOffset = sub_vecs[0].globalOffset();
  TEUCHOS_TEST_FOR_EXCEPT(globalOffset<0);
  const_iter_t v0_val = sub_vecs[0].values().begin();
  const ptrdiff_t v0_s = sub_vecs[0].stride();
  
  if( u_ < globalOffset || globalOffset + subDim - 1 < l_ ) {
    // None of the sub-vector elements that we are looking for is in this
    // vector chunk!
    return;
  }
 
  index_type
    i_l = ( l_ <= globalOffset ? 0 : l_ - globalOffset ),
    i_u = ( u_ >= globalOffset+subDim-1 ? subDim-1 : u_ - globalOffset );

  DRTSVV &drtsvv_reduct_obj = dyn_cast<DRTSVV>(*reduct_obj);
  SubVectorView<Scalar> sub_vec_targ = drtsvv_reduct_obj.get();

  const ArrayRCP<Scalar> svt_values = sub_vec_targ.values();

  for( index_type i = i_l; i <= i_u; ++i ) {
    svt_values[i+(globalOffset-l_)] = v0_val[i*v0_s];
  }

  drtsvv_reduct_obj.set(sub_vec_targ);

}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_GET_SUB_VECTOR_DEF_HPP
