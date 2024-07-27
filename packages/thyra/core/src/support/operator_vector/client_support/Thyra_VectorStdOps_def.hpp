// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_VECTOR_STD_OPS_HPP
#define THYRA_VECTOR_STD_OPS_HPP

#include "Thyra_VectorStdOps_decl.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "RTOpPack_ROpGetElement.hpp"
#include "RTOpPack_TOpSetElement.hpp"
#include "RTOpPack_ROpMin.hpp"
#include "RTOpPack_ROpMinIndex.hpp"
#include "RTOpPack_ROpMinIndexGreaterThanBound.hpp"
#include "RTOpPack_ROpMax.hpp"
#include "RTOpPack_ROpMaxIndex.hpp"
#include "RTOpPack_ROpMaxIndexLessThanBound.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_TOpAddScalar.hpp"
#include "RTOpPack_TOpEleWiseDivide.hpp"
#include "RTOpPack_TOpEleWiseProd.hpp"
#include "RTOpPack_TOpPairWiseMax.hpp"
#include "RTOpPack_TOpEleWiseConjProd.hpp"
#include "RTOpPack_TOpEleWiseProdUpdate.hpp"
#include "RTOpPack_TOpPairWiseMaxUpdate.hpp"
#include "RTOpPack_TOpRandomize.hpp"
#include "Teuchos_Assert.hpp"
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
  return v_rhs.norm_1();
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_2( const VectorBase<Scalar>& v_rhs )
{
  return v_rhs.norm_2();
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_2( const VectorBase<Scalar>& w, const VectorBase<Scalar>& v )
{
  return v.norm_2(w);
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_inf( const VectorBase<Scalar>& v_rhs )
{
  return v_rhs.norm_inf();
}


template<class Scalar>
Scalar Thyra::dot( const VectorBase<Scalar>& v_rhs1, const VectorBase<Scalar>& v_rhs2 )
{
  return v_rhs2.dot(v_rhs1);
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
  v_lhs->assign(alpha);
}


template<class Scalar>
void Thyra::copy( const VectorBase<Scalar>& v_rhs,
  const Ptr<VectorBase<Scalar> > &v_lhs )
{
  v_lhs->assign(v_rhs);
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
  v_lhs->scale(alpha);
}


template<class Scalar>
void Thyra::abs( const VectorBase<Scalar>& x, const Ptr<VectorBase<Scalar> > &y )
{
  y->abs(x);
}


template<class Scalar>
void Thyra::reciprocal( const VectorBase<Scalar>& x, const Ptr<VectorBase<Scalar> > &y )
{
  y->reciprocal(x);
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
void Thyra::pair_wise_max(
  const Scalar &alpha, const VectorBase<Scalar>& v_rhs1,
  const VectorBase<Scalar>& v_rhs2, const Ptr<VectorBase<Scalar> > &v_lhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpPairWiseMax<Scalar> pair_wise_max_op(alpha);
  applyOp<Scalar>( pair_wise_max_op, tuple(ptrInArg(v_rhs1),ptrInArg(v_rhs2)),
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
void Thyra::ele_wise_scale( const VectorBase<Scalar>& x,
  const Ptr<VectorBase<Scalar> > &y )
{
  y->ele_wise_scale(x);
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
void Thyra::pair_wise_max_update(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1,
  const Ptr<VectorBase<Scalar> > &v_lhs
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpPairWiseMaxUpdate<Scalar> pair_wise_max_update_op(alpha);
  applyOp<Scalar>( pair_wise_max_update_op, tuple(ptrInArg(v_rhs1)),
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
  y->linear_combination(alpha, x, beta);
}


template<class Scalar>
void Thyra::seed_randomize( unsigned int s )
{
  RTOpPack::TOpRandomize<Scalar>::set_static_seed(s);
}


template<class Scalar>
void Thyra::randomize( Scalar l, Scalar u, const Ptr<VectorBase<Scalar> > &v )
{
  v->randomize(l, u);
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
  v_lhs->update(alpha, v_rhs);
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
  template void abs( const VectorBase< SCALAR > &x, const Ptr<VectorBase< SCALAR > > &y );  \
   \
  template void reciprocal( const VectorBase< SCALAR > &x, const Ptr<VectorBase< SCALAR > > &y );  \
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
  template void ele_wise_scale( const VectorBase<SCALAR>& x, \
    const Ptr<VectorBase<SCALAR> > &y ); \
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
  template void pair_wise_max(  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const VectorBase<SCALAR >& v_rhs2, const Ptr<VectorBase<SCALAR > > &v_lhs  \
    );  \
   \
  template void pair_wise_max_update(  \
    const SCALAR& alpha, const VectorBase<SCALAR >& v_rhs1,  \
    const Ptr<VectorBase<SCALAR > > &v_lhs  \
    );  \
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
