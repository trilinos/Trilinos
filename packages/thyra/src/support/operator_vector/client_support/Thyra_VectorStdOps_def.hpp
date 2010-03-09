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

#ifndef THYRA_VECTOR_STD_OPS_HPP
#define THYRA_VECTOR_STD_OPS_HPP

#include "Thyra_VectorStdOps_decl.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_ROpGetElement.hpp"
#include "RTOpPack_TOpSetElement.hpp"
#include "RTOpPack_ROpMin.hpp"
#include "RTOpPack_ROpMinIndex.hpp"
#include "RTOpPack_ROpMinIndexGreaterThanBound.hpp"
#include "RTOpPack_ROpMax.hpp"
#include "RTOpPack_ROpMaxIndex.hpp"
#include "RTOpPack_ROpMaxIndexLessThanBound.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNorm2.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_ROpWeightedNorm2.hpp"
#include "RTOpPack_TOpAbs.hpp"
#include "RTOpPack_TOpAddScalar.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
#include "RTOpPack_TOpAssignVectors.hpp"
#include "RTOpPack_TOpAXPY.hpp"
#include "RTOpPack_TOpEleWiseDivide.hpp"
#include "RTOpPack_TOpEleWiseProd.hpp"
#include "RTOpPack_TOpEleWiseConjProd.hpp"
#include "RTOpPack_TOpEleWiseProdUpdate.hpp"
#include "RTOpPack_TOpLinearCombination.hpp"
#include "RTOpPack_TOpScaleVector.hpp"
#include "RTOpPack_TOpReciprocal.hpp"
#include "RTOpPack_TOpRandomize.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"


//
// All scalar types
//


// Standard text names


// Reduction operations


template<class Scalar>
Scalar Thyra::sum( const VectorBase<Scalar>& v_rhs )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RCP<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
  applyOp<Scalar>(sum_op,
    tuple(ptrInArg(v_rhs)),
    ArrayView<Ptr<VectorBase<Scalar> > >(null),
    sum_targ.ptr() );
  return sum_op(*sum_targ);
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_1( const VectorBase<Scalar>& v_rhs )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpNorm1<Scalar> norm_1_op;
  Teuchos::RCP<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
  applyOp<Scalar>(norm_1_op, tuple(ptrInArg(v_rhs)),
    ArrayView<Ptr<VectorBase<Scalar> > >(null),
    norm_1_targ.ptr() );
  return norm_1_op(*norm_1_targ);
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_2( const VectorBase<Scalar>& v_rhs )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpNorm2<Scalar> norm_2_op;
  Teuchos::RCP<RTOpPack::ReductTarget> norm_2_targ = norm_2_op.reduct_obj_create();
  applyOp<Scalar>(norm_2_op, tuple(ptrInArg(v_rhs)),
    ArrayView<Ptr<VectorBase<Scalar> > >(null),
    norm_2_targ.ptr() );
  return norm_2_op(*norm_2_targ);
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_2( const VectorBase<Scalar>& w, const VectorBase<Scalar>& v )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpWeightedNorm2<Scalar> wght_norm_2_op;
  Teuchos::RCP<RTOpPack::ReductTarget> wght_norm_2_targ = wght_norm_2_op.reduct_obj_create();
  applyOp<Scalar>(wght_norm_2_op, tuple(ptrInArg(w), ptrInArg(v)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    wght_norm_2_targ.ptr());
  return wght_norm_2_op(*wght_norm_2_targ);
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_inf( const VectorBase<Scalar>& v_rhs )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpNormInf<Scalar> norm_inf_op;
  Teuchos::RCP<RTOpPack::ReductTarget> norm_inf_targ = norm_inf_op.reduct_obj_create();
  applyOp<Scalar>(norm_inf_op, tuple(ptrInArg(v_rhs)),
    ArrayView<Ptr<VectorBase<Scalar> > >(null),
    norm_inf_targ.ptr() );
  return norm_inf_op(*norm_inf_targ);
}


template<class Scalar>
Scalar Thyra::dot( const VectorBase<Scalar>& v_rhs1, const VectorBase<Scalar>& v_rhs2 )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpDotProd<Scalar> dot_prod_op;
  Teuchos::RCP<RTOpPack::ReductTarget>
    dot_prod_targ = dot_prod_op.reduct_obj_create();
  applyOp<Scalar>(dot_prod_op,
    tuple(ptrInArg(v_rhs1), ptrInArg(v_rhs2))(),
    ArrayView<Ptr<VectorBase<Scalar> > >(null),
    dot_prod_targ.ptr()
    );
  return dot_prod_op(*dot_prod_targ);
}


template<class Scalar>
Scalar Thyra::get_ele( const VectorBase<Scalar>& v, Ordinal i )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(i, 0, v.space()->dim());
#endif
  RTOpPack::ROpGetElement<Scalar> get_ele_op(i);
  Teuchos::RCP<RTOpPack::ReductTarget> get_ele_targ = get_ele_op.reduct_obj_create();
  applyOp<Scalar>(get_ele_op, tuple(ptrInArg(v)),
    ArrayView<Ptr<VectorBase<Scalar> > >(null),
    get_ele_targ.ptr() );
  return get_ele_op(*get_ele_targ);
}


// Transformation operations


template<class Scalar>
void Thyra::set_ele( Ordinal i, Scalar alpha, const Ptr<VectorBase<Scalar> > &v )
{
  using Teuchos::tuple; using Teuchos::null;
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(i, 0, v->space()->dim());
#endif
  RTOpPack::TOpSetElement<Scalar> set_ele_op(i, alpha);
  applyOp<Scalar>(set_ele_op,
    ArrayView<Ptr<const VectorBase<Scalar> > >(null),
    tuple(v),
    null);
}


template<class Scalar>
void Thyra::put_scalar( const Scalar& alpha, const Ptr<VectorBase<Scalar> > &v_lhs )
{
  using Teuchos::tuple; using Teuchos::null;
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
  applyOp<Scalar>(assign_scalar_op,
    ArrayView<Ptr<const VectorBase<Scalar> > >(null),
    tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::copy( const VectorBase<Scalar>& v_rhs,
  const Ptr<VectorBase<Scalar> > &v_lhs )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpAssignVectors<Scalar> assign_vectors_op;
  applyOp<Scalar>( assign_vectors_op, tuple(ptrInArg(v_rhs)), tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::add_scalar( const Scalar& alpha, const Ptr<VectorBase<Scalar> > &v_lhs )
{
  using Teuchos::tuple; using Teuchos::null;
  RTOpPack::TOpAddScalar<Scalar> add_scalar_op(alpha);
  applyOp<Scalar>(add_scalar_op,
    ArrayView<Ptr<const VectorBase<Scalar> > >(null),
    tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::scale( const Scalar& alpha, const Ptr<VectorBase<Scalar> > &v_lhs )
{
  using Teuchos::tuple; using Teuchos::null;
  if( alpha == ScalarTraits<Scalar>::zero() ) {
    assign(v_lhs, ScalarTraits<Scalar>::zero());
  }
  else if( alpha != ScalarTraits<Scalar>::one() ) {
    RTOpPack::TOpScaleVector<Scalar> scale_vector_op(alpha);
    applyOp<Scalar>(scale_vector_op,
      ArrayView<Ptr<const VectorBase<Scalar> > >(null),
      tuple(v_lhs), null );
  }
}


template<class Scalar>
void Thyra::abs( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpAbs<Scalar> abs_op;
  applyOp<Scalar>( abs_op, tuple(ptrInArg(x)), tuple(y), null );
}


template<class Scalar>
void Thyra::reciprocal( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpReciprocal<Scalar> recip_op;
  applyOp<Scalar>( recip_op, tuple(ptrInArg(x)), tuple(y), null );
}


template<class Scalar>
void Thyra::ele_wise_prod(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1,
  const VectorBase<Scalar>& v_rhs2, const Ptr<VectorBase<Scalar> > &v_lhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpEleWiseProd<Scalar> ele_wise_prod_op(alpha);
  applyOp<Scalar>( ele_wise_prod_op, tuple(ptrInArg(v_rhs1),ptrInArg(v_rhs2)),
    tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::ele_wise_conj_prod(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1,
  const VectorBase<Scalar>& v_rhs2, const Ptr<VectorBase<Scalar> > &v_lhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpEleWiseConjProd<Scalar> ele_wise_conj_prod_op(alpha);
  applyOp<Scalar>( ele_wise_conj_prod_op, tuple(ptrInArg(v_rhs1),ptrInArg(v_rhs2)),
    tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::Vp_StVtV(
  const Ptr<VectorBase<Scalar> > &v_lhs,
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1,
  const VectorBase<Scalar>& v_rhs2
  )
{
  ele_wise_prod(alpha,v_rhs1,v_rhs2,v_lhs);
}


template<class Scalar>
void Thyra::ele_wise_prod_update(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1,
  const Ptr<VectorBase<Scalar> > &v_lhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpEleWiseProdUpdate<Scalar> ele_wise_prod_update_op(alpha);
  applyOp<Scalar>( ele_wise_prod_update_op, tuple(ptrInArg(v_rhs1)),
    tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::Vt_StV(
  const Ptr<VectorBase<Scalar> > &v_lhs,
  const Scalar& alpha, const VectorBase<Scalar>& x )
{
  ele_wise_prod_update(alpha,x,v_lhs);
}


template<class Scalar>
void Thyra::ele_wise_divide(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1,
  const VectorBase<Scalar>& v_rhs2,
  const Ptr<VectorBase<Scalar> > &v_lhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpEleWiseDivide<Scalar> ele_wise_divide_op(alpha);
  applyOp<Scalar>( ele_wise_divide_op, tuple(ptrInArg(v_rhs1),ptrInArg(v_rhs2)),
    tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::linear_combination(
  const ArrayView<const Scalar> &alpha,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &x,
  const Scalar &beta,
  const Ptr<VectorBase<Scalar> > &y
  )
{
  using Teuchos::tuple; using Teuchos::ptr; using Teuchos::ptrInArg;
  using Teuchos::null;
  const int m = x.size();
  if( beta == Teuchos::ScalarTraits<Scalar>::one() && m == 1 ) {
    Vp_StV( y, alpha[0], *x[0] );
    return;
  }
  else if( m == 0 ) {
    Vt_S( y, beta );
    return;
  }
  RTOpPack::TOpLinearCombination<Scalar> lin_comb_op(alpha,beta);
  applyOp<Scalar>( lin_comb_op, x, tuple(y), null );
}


template<class Scalar>
void Thyra::seed_randomize( unsigned int s )
{
  RTOpPack::TOpRandomize<Scalar>::set_static_seed(s);
}


template<class Scalar>
void Thyra::randomize( Scalar l, Scalar u, const Ptr<VectorBase<Scalar> > &v )
{
  using Teuchos::tuple; using Teuchos::null;
  RTOpPack::TOpRandomize<Scalar> random_vector_op(l,u);
  applyOp<Scalar>( random_vector_op,
    ArrayView<Ptr<const VectorBase<Scalar> > >(null),
    tuple(v),
    null );
  // Warning! If the RTOpPack::TOpRandomize<Scalar> object is ever made
  // static, the one must be careful to change the seed in between calls.
  // Right now the seed is being incremented by the constructor automatically.
  // It is important to generate different random vectors on each call
  // (i.e. to generate different columns in a multi-vector).
}


// Linear algebra names


template<class Scalar>
void Thyra::assign( const Ptr<VectorBase<Scalar> > &v_lhs, const Scalar& alpha )
{
  put_scalar(alpha,v_lhs);
}


template<class Scalar>
void Thyra::assign( const Ptr<VectorBase<Scalar> > &v_lhs, const VectorBase<Scalar>& v_rhs )
{
  copy(v_rhs,v_lhs);
}


template<class Scalar>
void Thyra::Vp_S( const Ptr<VectorBase<Scalar> > &v_lhs, const Scalar& alpha )
{
  add_scalar(alpha,v_lhs);
}


template<class Scalar>
void Thyra::Vt_S( const Ptr<VectorBase<Scalar> > &v_lhs, const Scalar& alpha )
{
  scale(alpha,v_lhs);
}


template<class Scalar>
void Thyra::V_StV( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha,
  const VectorBase<Scalar> &x
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  linear_combination<Scalar>( tuple<Scalar>(alpha), tuple(ptrInArg(x)),
    ScalarTraits<Scalar>::zero(), y );
}


template<class Scalar>
void Thyra::Vp_StV( const Ptr<VectorBase<Scalar> > &v_lhs, const Scalar& alpha,
  const VectorBase<Scalar>& v_rhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpAXPY<Scalar> axpy_op(alpha);
  applyOp<Scalar>( axpy_op, tuple(ptrInArg(v_rhs)), tuple(v_lhs), null );
}


template<class Scalar>
void Thyra::Vp_V( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x,
  const Scalar& beta
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  linear_combination<Scalar>(
    tuple<Scalar>(Teuchos::ScalarTraits<Scalar>::one()),
    tuple(ptrInArg(x)),
    beta, y );
}


template<class Scalar>
void Thyra::V_V( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x )
{
  assign(y,x);
}


template<class Scalar>
void Thyra::V_S( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha )
{
  assign(y,alpha);
}


template<class Scalar>
void Thyra::V_VpV( const Ptr<VectorBase<Scalar> > &z, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& y
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(ST::one(),ST::one()),
    tuple(ptrInArg(x),ptrInArg(y)),
    ST::zero(), z
    );
}


template<class Scalar>
void Thyra::V_VmV( const Ptr<VectorBase<Scalar> > &z, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& y
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(ST::one(),Scalar(-ST::one())),
    tuple(ptrInArg(x),ptrInArg(y)),
    ST::zero(), z
    );
}


template<class Scalar>
void Thyra::V_StVpV( const Ptr<VectorBase<Scalar> > &z, const Scalar &alpha,
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& y
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(alpha, ST::one()), tuple(ptrInArg(x),ptrInArg(y)),
    ST::zero(), z
    );
}


template<class Scalar>
void Thyra::V_VpStV( const Ptr<VectorBase<Scalar> > &z,
  const VectorBase<Scalar>& x,
  const Scalar &alpha, const VectorBase<Scalar>& y )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(ST::one(), alpha), tuple(ptrInArg(x),ptrInArg(y)),
    ST::zero(), z
    );
}


template<class Scalar>
void Thyra::V_StVpStV( const Ptr<VectorBase<Scalar> > &z, const Scalar &alpha,
  const VectorBase<Scalar>& x, const Scalar &beta, const VectorBase<Scalar>& y
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(alpha, beta), tuple(ptrInArg(x),ptrInArg(y)),
    ST::zero(), z
    );
}


//
// For real types only
//


template<class Scalar>
Scalar Thyra::min( const VectorBase<Scalar>& x ) {
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpMin<Scalar> min_op;
  Teuchos::RCP<RTOpPack::ReductTarget> min_targ = min_op.reduct_obj_create();
  applyOp<Scalar>( min_op, tuple(ptrInArg(x)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    min_targ.ptr() );
  return min_op(*min_targ);
}


template<class Scalar>
void Thyra::min( const VectorBase<Scalar>& x,
  const Ptr<Scalar> &minEle, const Ptr<Ordinal> &minIndex
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpMinIndex<Scalar> min_op;
  Teuchos::RCP<RTOpPack::ReductTarget> min_targ = min_op.reduct_obj_create();
  applyOp<Scalar>( min_op, tuple(ptrInArg(x)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    min_targ.ptr() );
  RTOpPack::ScalarIndex<Scalar> scalarIndex = min_op(*min_targ);
  *minEle = scalarIndex.scalar;
  *minIndex = scalarIndex.index;
}


template<class Scalar>
void Thyra::minGreaterThanBound( const VectorBase<Scalar>& x,
  const Scalar &bound, const Ptr<Scalar> &minEle, const Ptr<Ordinal> &minIndex
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpMinIndexGreaterThanBound<Scalar> min_op(bound);
  Teuchos::RCP<RTOpPack::ReductTarget> min_targ = min_op.reduct_obj_create();
  applyOp<Scalar>( min_op, tuple(ptrInArg(x)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    min_targ.ptr() );
  RTOpPack::ScalarIndex<Scalar> scalarIndex = min_op(*min_targ);
  *minEle = scalarIndex.scalar;
  *minIndex = scalarIndex.index;
}


template<class Scalar>
Scalar Thyra::max( const VectorBase<Scalar>& x )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpMax<Scalar> max_op;
  Teuchos::RCP<RTOpPack::ReductTarget> max_targ = max_op.reduct_obj_create();
  applyOp<Scalar>( max_op, tuple(ptrInArg(x)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    max_targ.ptr() );
  return max_op(*max_targ);
}


template<class Scalar>
void Thyra::max( const VectorBase<Scalar>& x,
  const Ptr<Scalar> &maxEle, const Ptr<Ordinal> &maxIndex
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpMaxIndex<Scalar> max_op;
  Teuchos::RCP<RTOpPack::ReductTarget> max_targ = max_op.reduct_obj_create();
  applyOp<Scalar>( max_op, tuple(ptrInArg(x)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    max_targ.ptr() );
  RTOpPack::ScalarIndex<Scalar> scalarIndex = max_op(*max_targ);
  *maxEle = scalarIndex.scalar;
  *maxIndex = scalarIndex.index;
}


template<class Scalar>
void Thyra::maxLessThanBound( const VectorBase<Scalar>& x,
  const Scalar &bound, const Ptr<Scalar> &maxEle, const Ptr<Ordinal> &maxIndex
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::ROpMaxIndexLessThanBound<Scalar> max_op(bound);
  Teuchos::RCP<RTOpPack::ReductTarget> max_targ = max_op.reduct_obj_create();
  applyOp<Scalar>( max_op, tuple(ptrInArg(x)),
    ArrayView<const Ptr<VectorBase<Scalar> > >(null),
    max_targ.ptr() );
  RTOpPack::ScalarIndex<Scalar> scalarIndex = max_op(*max_targ);
  *maxEle = scalarIndex.scalar;
  *maxIndex = scalarIndex.index;
}


//
// Explicit instantiation macro
//


#define THYRA_VECTOR_STD_OPS_INSTANT(SCALAR) \
   \
  template SCALAR sum( const VectorBase<SCALAR >& v_rhs );  \
   \
  template ScalarTraits<SCALAR >::magnitudeType  \
  norm_1( const VectorBase<SCALAR >& v_rhs );  \
   \
  template ScalarTraits<SCALAR >::magnitudeType  \
  norm_2( const VectorBase<SCALAR >& v_rhs );  \
   \
  template ScalarTraits<SCALAR >::magnitudeType  \
  norm_2( const VectorBase<SCALAR >& w, const VectorBase<SCALAR >& v );  \
   \
  template ScalarTraits<SCALAR >::magnitudeType  \
  norm_inf( const VectorBase<SCALAR >& v_rhs );  \
   \
  template SCALAR dot( const VectorBase<SCALAR >& v_rhs1, const VectorBase<SCALAR >& v_rhs2 );  \
   \
  template SCALAR get_ele( const VectorBase<SCALAR >& v, Ordinal i );  \
   \
  template void set_ele( Ordinal i, SCALAR alpha, const Ptr<VectorBase<SCALAR > > &v );  \
   \
  template void put_scalar( const SCALAR& alpha, const Ptr<VectorBase<SCALAR > > &v_lhs );  \
   \
  template void copy( const VectorBase<SCALAR >& v_rhs,  \
    const Ptr<VectorBase<SCALAR > > &v_lhs );  \
   \
  template void add_scalar( const SCALAR& alpha, const Ptr<VectorBase<SCALAR > > &v_lhs );  \
   \
  template void scale( const SCALAR& alpha, const Ptr<VectorBase<SCALAR > > &v_lhs );  \
   \
  template void abs( const Ptr<VectorBase<SCALAR > > &y, const VectorBase<SCALAR >& x );  \
   \
  template void reciprocal( const Ptr<VectorBase<SCALAR > > &y, const VectorBase<SCALAR >& x );  \
   \
  template void ele_wise_prod(  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const VectorBase<SCALAR >& v_rhs2, const Ptr<VectorBase<SCALAR > > &v_lhs  \
    );  \
   \
  template void ele_wise_conj_prod(  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const VectorBase<SCALAR >& v_rhs2, const Ptr<VectorBase<SCALAR > > &v_lhs  \
    );  \
   \
  template void Vp_StVtV(  \
    const Ptr<VectorBase<SCALAR > > &v_lhs,  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const VectorBase<SCALAR >& v_rhs2  \
    );  \
   \
  template void ele_wise_prod_update(  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const Ptr<VectorBase<SCALAR > > &v_lhs  \
    );  \
   \
  template void Vt_StV(  \
    const Ptr<VectorBase<SCALAR > > &v_lhs,  \
    const SCALAR& alpha, const VectorBase<SCALAR >& x );  \
   \
  template void ele_wise_divide(  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const VectorBase<SCALAR >& v_rhs2,  \
    const Ptr<VectorBase<SCALAR > > &v_lhs  \
    );  \
   \
  template void linear_combination(  \
    const ArrayView<const SCALAR > &alpha,  \
    const ArrayView<const Ptr<const VectorBase<SCALAR > > > &x,  \
    const SCALAR &beta,  \
    const Ptr<VectorBase<SCALAR > > &y  \
    );  \
   \
  template void seed_randomize<SCALAR >( unsigned int s );  \
   \
  template void randomize( SCALAR l, SCALAR u, const Ptr<VectorBase<SCALAR > > &v );  \
   \
  template void assign( const Ptr<VectorBase<SCALAR > > &v_lhs, const SCALAR& alpha );  \
   \
  template void assign( const Ptr<VectorBase<SCALAR > > &v_lhs, const VectorBase<SCALAR >& v_rhs );  \
   \
  template void Vp_S( const Ptr<VectorBase<SCALAR > > &v_lhs, const SCALAR& alpha );  \
   \
  template void Vt_S( const Ptr<VectorBase<SCALAR > > &v_lhs, const SCALAR& alpha );  \
   \
  template void V_StV( const Ptr<VectorBase<SCALAR > > &y, const SCALAR& alpha,  \
    const VectorBase<SCALAR > &x  \
    );  \
   \
  template void Vp_StV( const Ptr<VectorBase<SCALAR > > &v_lhs, const SCALAR& alpha,  \
    const VectorBase<SCALAR >& v_rhs  \
    );  \
   \
  template void Vp_V( const Ptr<VectorBase<SCALAR > > &y, const VectorBase<SCALAR >& x,  \
    const SCALAR& beta  \
    );  \
   \
  template void V_V( const Ptr<VectorBase<SCALAR > > &y, const VectorBase<SCALAR >& x );  \
   \
  template void V_S( const Ptr<VectorBase<SCALAR > > &y, const SCALAR& alpha );  \
   \
  template void V_VpV( const Ptr<VectorBase<SCALAR > > &z, const VectorBase<SCALAR >& x,  \
    const VectorBase<SCALAR >& y  \
    );  \
   \
  template void V_VmV( const Ptr<VectorBase<SCALAR > > &z, const VectorBase<SCALAR >& x,  \
    const VectorBase<SCALAR >& y  \
    );  \
   \
  template void V_StVpV( const Ptr<VectorBase<SCALAR > > &z, const SCALAR &alpha,  \
    const VectorBase<SCALAR >& x, const VectorBase<SCALAR >& y  \
    );  \
   \
  template void V_VpStV( const Ptr<VectorBase<SCALAR > > &z,  \
    const VectorBase<SCALAR >& x,  \
    const SCALAR &alpha, const VectorBase<SCALAR >& y );  \
   \
  template void V_StVpStV( const Ptr<VectorBase<SCALAR > > &z, const SCALAR &alpha,  \
    const VectorBase<SCALAR >& x, const SCALAR &beta, const VectorBase<SCALAR >& y  \
    );  \



#define THYRA_VECTOR_STD_OPS_REAL_INSTANT(SCALAR) \
   \
  template SCALAR min( const VectorBase<SCALAR >& x );  \
   \
  template void min( const VectorBase<SCALAR >& x,  \
    const Ptr<SCALAR > &minEle, const Ptr<Ordinal> &minIndex  \
    );  \
   \
  template void minGreaterThanBound( const VectorBase<SCALAR >& x,  \
    const SCALAR &bound, const Ptr<SCALAR > &minEle, const Ptr<Ordinal> &minIndex  \
    );  \
   \
  template SCALAR max( const VectorBase<SCALAR >& x );  \
   \
  template void max( const VectorBase<SCALAR >& x,  \
    const Ptr<SCALAR > &maxEle, const Ptr<Ordinal> &maxIndex  \
    );  \
   \
  template void maxLessThanBound( const VectorBase<SCALAR >& x,  \
    const SCALAR &bound, const Ptr<SCALAR > &maxEle, const Ptr<Ordinal> &maxIndex  \
    );  \


#endif // THYRA_VECTOR_STD_OPS_HPP
