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
void ROpGetSubVector<Scalar>::get_reduct_type_num_entries(
  int* num_values
  ,int* num_indexes
  ,int* num_chars
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
ROpGetSubVector<Scalar>::reduct_obj_create() const
{
  const index_type subDim = u_ - l_ + 1;
  return defaultReductTarget(
    SubVectorView<Scalar>( l_, subDim, Teuchos::arcp<Scalar>(subDim), 1 )
    );
}


template<class Scalar>
void ROpGetSubVector<Scalar>::reduce_reduct_objs(
  const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
  ) const
{

  using Teuchos::dyn_cast;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;

  DRTSVV &drtsvv_inout_reduct_obj = dyn_cast<DRTSVV>(*inout_reduct_obj);

  const ConstSubVectorView<Scalar> sub_vec_in =
    dyn_cast<const DRTSVV>(in_reduct_obj).get();
  SubVectorView<Scalar> sub_vec_inout = drtsvv_inout_reduct_obj.get();

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(
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
void ROpGetSubVector<Scalar>::reduct_obj_reinit( ReductTarget* reduct_obj ) const
{
  using Teuchos::dyn_cast;
  typedef typename ArrayRCP<Scalar>::iterator iter_t;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;
  DRTSVV &drtsvv_inout_reduct_obj = dyn_cast<DRTSVV>(*reduct_obj);
  SubVectorView<Scalar> sub_vec = drtsvv_inout_reduct_obj.get();
  std::fill( sub_vec.values().begin(), sub_vec.values().end(),
    ScalarTraits<Scalar>::zero() );
}


template<class Scalar>
void ROpGetSubVector<Scalar>::extract_reduct_obj_state(
  const ReductTarget &reduct_obj
  ,int num_values
  ,primitive_value_type value_data[]
  ,int num_indexes
  ,Teuchos_Index index_data[]
  ,int num_chars
  ,::RTOpPack::char_type char_data[]
  ) const
{
  using Teuchos::null;
  using Teuchos::dyn_cast;
  using Teuchos::arrayView;
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
      arrayView(value_data+value_data_off, num_prim_objs_per_scalar ),
      null, null );
  }
}


template<class Scalar>
void ROpGetSubVector<Scalar>::load_reduct_obj_state(
  int num_values
  ,const primitive_value_type value_data[]
  ,int num_indexes
  ,const Teuchos_Index index_data[]
  ,int num_chars
  ,const ::RTOpPack::char_type char_data[]
  ,ReductTarget *reduct_obj
  ) const
{
  using Teuchos::null;
  using Teuchos::ptr;
  using Teuchos::outArg;
  using Teuchos::dyn_cast;
  using Teuchos::arrayView;
  using Teuchos::arcp_const_cast;
  typedef PrimitiveTypeTraits<Scalar,Scalar> PTT;
  typedef DefaultReductTarget<SubVectorView<Scalar> > DRTSVV;
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  DRTSVV &drtsvv_reduct_obj = dyn_cast<DRTSVV>(*ptr(reduct_obj));
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
      arrayView(value_data+value_data_off, num_prim_objs_per_scalar), null, null,
      outArg(sv_values[k]) );
  }
}


template<class Scalar>
void ROpGetSubVector<Scalar>::get_op_type_num_entries(
  int* num_values
  ,int* num_indexes
  ,int* num_chars
  ) const
{
  TEST_FOR_EXCEPT(true); // ToDo: Implement!
}


template<class Scalar>
void ROpGetSubVector<Scalar>::extract_op_state(
  int num_values
  ,primitive_value_type value_data[]
  ,int num_indexes
  ,Teuchos_Index index_data[]
  ,int num_chars
  ,::RTOpPack::char_type char_data[]
  ) const
{
  TEST_FOR_EXCEPT(true); // ToDo: Implement!
}


template<class Scalar>
void ROpGetSubVector<Scalar>::load_op_state(
  int num_values
  ,const primitive_value_type value_data[]
  ,int num_indexes
  ,const Teuchos_Index index_data[]
  ,int num_chars
  ,const ::RTOpPack::char_type char_data[]
  )
{
  TEST_FOR_EXCEPT(true); // ToDo: Implement!
}


template<class Scalar>
bool ROpGetSubVector<Scalar>::coord_invariant() const
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
    sub_vecs, targ_sub_vecs, reduct_obj );

  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
  const index_type subDim  = sub_vecs[0].subDim();
  const index_type globalOffset = sub_vecs[0].globalOffset();
  TEST_FOR_EXCEPT(globalOffset<0);
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
