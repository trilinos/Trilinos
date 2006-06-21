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

#include "Thyra_VectorStdOpsDecl.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "RTOpPack_ROpDotProd.hpp"
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
#include "RTOpPack_TOpEleWiseProdUpdate.hpp"
#include "RTOpPack_TOpLinearCombination.hpp"
#include "RTOpPack_TOpScaleVector.hpp"
#include "RTOpPack_TOpReciprocal.hpp"
#include "RTOpPack_TOpRandomize.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_arrayArg.hpp"

//
// All scalar types
//

// Standard test names

// Reduction operations

template<class Scalar>
Scalar Thyra::sum( const VectorBase<Scalar>& v_rhs )
{
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &v_rhs };
  applyOp<Scalar>(sum_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*sum_targ);
  return sum_op(*sum_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_1( const VectorBase<Scalar>& v_rhs )
{
  RTOpPack::ROpNorm1<Scalar> norm_1_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &v_rhs };
  applyOp<Scalar>(norm_1_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*norm_1_targ);
  return norm_1_op(*norm_1_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_2( const VectorBase<Scalar>& v_rhs )
{
  RTOpPack::ROpNorm2<Scalar> norm_2_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_2_targ = norm_2_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &v_rhs };
  applyOp<Scalar>(norm_2_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*norm_2_targ);
  return norm_2_op(*norm_2_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_2( const VectorBase<Scalar>& w, const VectorBase<Scalar>& v )
{
  RTOpPack::ROpWeightedNorm2<Scalar> wght_norm_2_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> wght_norm_2_targ = wght_norm_2_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &w, &v };
  applyOp<Scalar>(wght_norm_2_op,2,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*wght_norm_2_targ);
  return wght_norm_2_op(*wght_norm_2_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_inf( const VectorBase<Scalar>& v_rhs )
{
  RTOpPack::ROpNormInf<Scalar> norm_inf_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_inf_targ = norm_inf_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &v_rhs };
  applyOp<Scalar>(norm_inf_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*norm_inf_targ);
  return norm_inf_op(*norm_inf_targ);
}

template<class Scalar>
Scalar Thyra::dot( const VectorBase<Scalar>& v_rhs1, const VectorBase<Scalar>& v_rhs2 )
{
  RTOpPack::ROpDotProd<Scalar> dot_prod_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> dot_prod_targ = dot_prod_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &v_rhs1, &v_rhs2 };
  applyOp<Scalar>(dot_prod_op,2,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*dot_prod_targ);
  return dot_prod_op(*dot_prod_targ);
}

template<class Scalar>
Scalar Thyra::get_ele( const VectorBase<Scalar>& v, Index i )
{
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &v };
  applyOp<Scalar>(sum_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*sum_targ,i,1,0);
  return sum_op(*sum_targ);
}

// Transformation operations

template<class Scalar>
void Thyra::set_ele( Index i, Scalar alpha, VectorBase<Scalar>* v )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v==NULL,std::logic_error,"set_ele(...), Error!");
#endif
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
  VectorBase<Scalar>* targ_vecs[] = { v };
  applyOp<Scalar>(assign_scalar_op,0,(const VectorBase<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL,i,1,0);
}

template<class Scalar>
void Thyra::put_scalar( const Scalar& alpha, VectorBase<Scalar>* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"put_scalar(...), Error!");
#endif
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
  VectorBase<Scalar>* targ_vecs[] = { v_lhs };
  applyOp<Scalar>(assign_scalar_op,0,(const VectorBase<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::copy( const VectorBase<Scalar>& v_rhs, VectorBase<Scalar>* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"copy(...), Error!");
#endif
  RTOpPack::TOpAssignVectors<Scalar> assign_vectors_op;
  const VectorBase<Scalar>* vecs[]      = { &v_rhs };
  VectorBase<Scalar>*       targ_vecs[] = { v_lhs  };
  applyOp<Scalar>(assign_vectors_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::add_scalar( const Scalar& alpha, VectorBase<Scalar>* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"add_scalar(...), Error!");
#endif
  RTOpPack::TOpAddScalar<Scalar> add_scalar_op(alpha);
  VectorBase<Scalar>* targ_vecs[] = { v_lhs };
  applyOp<Scalar>(add_scalar_op,0,(const VectorBase<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::scale( const Scalar& alpha, VectorBase<Scalar>* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"scale(...), Error!");
#endif
  if( alpha == Teuchos::ScalarTraits<Scalar>::zero() ) {
    assign(v_lhs,Teuchos::ScalarTraits<Scalar>::zero());
  }
  else if( alpha != Teuchos::ScalarTraits<Scalar>::one() ) {
    RTOpPack::TOpScaleVector<Scalar> scale_vector_op(alpha);
    VectorBase<Scalar>* targ_vecs[] = { v_lhs };
    applyOp<Scalar>(scale_vector_op,0,(const VectorBase<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  }
}

template<class Scalar>
void Thyra::abs( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(y==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAbs<Scalar> abs_op;
  const VectorBase<Scalar>* vecs[]      = { &x };
  VectorBase<Scalar>*       targ_vecs[] = { y  };
  applyOp<Scalar>(abs_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::reciprocal( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(y==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpReciprocal<Scalar> recip_op;
  const VectorBase<Scalar>* vecs[]      = { &x };
  VectorBase<Scalar>*       targ_vecs[] = { y  };
  applyOp<Scalar>(recip_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::ele_wise_prod(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1, const VectorBase<Scalar>& v_rhs2
  ,VectorBase<Scalar>* v_lhs
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_prod(...), Error");
#endif
  RTOpPack::TOpEleWiseProd<Scalar> ele_wise_prod_op(alpha);
  const VectorBase<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
  VectorBase<Scalar>*       targ_vecs[] = { v_lhs };
  applyOp<Scalar>(ele_wise_prod_op,2,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::Vp_StVtV(
  VectorBase<Scalar>* v_lhs
  ,const Scalar& alpha, const VectorBase<Scalar>& v_rhs1, const VectorBase<Scalar>& v_rhs2
  )
{
  ele_wise_prod(alpha,v_rhs1,v_rhs2,v_lhs);
}

template<class Scalar>
void Thyra::ele_wise_prod_update(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1, VectorBase<Scalar>* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_prod_update(...), Error");
#endif
  RTOpPack::TOpEleWiseProdUpdate<Scalar> ele_wise_prod_update_op(alpha);
  const VectorBase<Scalar>* vecs[]      = { &v_rhs1 };
  VectorBase<Scalar>*       targ_vecs[] = { v_lhs };
  applyOp<Scalar>(ele_wise_prod_update_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::Vt_StV(
  VectorBase<Scalar>* v_lhs, const Scalar& alpha, const VectorBase<Scalar>& x )
{
  ele_wise_prod_update(alpha,x,v_lhs);
}

template<class Scalar>
void Thyra::ele_wise_divide(
  const Scalar& alpha, const VectorBase<Scalar>& v_rhs1, const VectorBase<Scalar>& v_rhs2
  ,VectorBase<Scalar>* v_lhs
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_divide(...), Error");
#endif
  RTOpPack::TOpEleWiseDivide<Scalar> ele_wise_divide_op(alpha);
  const VectorBase<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
  VectorBase<Scalar>*       targ_vecs[] = { v_lhs };
  applyOp<Scalar>(ele_wise_divide_op,2,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::linear_combination(
  const int                    m
  ,const Scalar                alpha[]
  ,const VectorBase<Scalar>*   x[]
  ,const Scalar                &beta
  ,VectorBase<Scalar>          *y
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(y==NULL,std::logic_error,"linear_combination(...), Error!");
#endif
  if( beta == Teuchos::ScalarTraits<Scalar>::one() && m == 1 ) {
    Vp_StV( y, alpha[0], *x[0] );
    return;
  }
  else if( m == 0 ) {
    Vt_S( y, beta );
    return;
  }
  RTOpPack::TOpLinearCombination<Scalar> lin_comb_op(m,alpha,beta);
  VectorBase<Scalar>* targ_vecs[] = { y };
  applyOp<Scalar>(lin_comb_op,m,x,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::seed_randomize( unsigned int s )
{
  RTOpPack::TOpRandomize<Scalar>::set_static_seed(s);
}

template<class Scalar>
void Thyra::randomize( Scalar l, Scalar u, VectorBase<Scalar>* v )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v==NULL,std::logic_error,"Vt_S(...), Error");
#endif
  RTOpPack::TOpRandomize<Scalar> random_vector_op(l,u);
  VectorBase<Scalar>* targ_vecs[] = { v };
  applyOp<Scalar>(random_vector_op,0,(const VectorBase<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  // Warning! If the RTOpPack::TOpRandomize<Scalar> object is ever made
  // static, the one must be careful to change the seed in between calls.
  // Right now the seed is being incremented by the constructor automatically.
  // It is important to generate different random vectors on each call
  // (i.e. to generate different columns in a multi-vector).
}

// Linear algebra names

template<class Scalar>
void Thyra::assign( VectorBase<Scalar>* v_lhs, const Scalar& alpha )
{
  put_scalar(alpha,v_lhs);
}

template<class Scalar>
void Thyra::assign( VectorBase<Scalar>* v_lhs, const VectorBase<Scalar>& v_rhs )
{
  copy(v_rhs,v_lhs);
}

template<class Scalar>
void Thyra::Vp_S( VectorBase<Scalar>* v_lhs, const Scalar& alpha )
{
  add_scalar(alpha,v_lhs);
}

template<class Scalar>
void Thyra::Vt_S(
  VectorBase<Scalar>* v_lhs, const Scalar& alpha )
{
  scale(alpha,v_lhs);
}

template<class Scalar>
void Thyra::V_StV( VectorBase<Scalar>* y, const Scalar& alpha, const VectorBase<Scalar> &x )
{
  linear_combination(
    1,Teuchos::arrayArg<Scalar>(alpha)(),Teuchos::arrayArg<const VectorBase<Scalar>*>(&x)()
    ,Teuchos::ScalarTraits<Scalar>::zero(),y
    );
}

template<class Scalar>
void Thyra::Vp_StV( VectorBase<Scalar>* v_lhs, const Scalar& alpha, const VectorBase<Scalar>& v_rhs )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vp_StV(...), Error!");
#endif
  RTOpPack::TOpAXPY<Scalar> axpy_op(alpha);
  const VectorBase<Scalar>* vecs[]      = { &v_rhs };
  VectorBase<Scalar>*       targ_vecs[] = { v_lhs  };
  applyOp<Scalar>(axpy_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void Thyra::Vp_V( VectorBase<Scalar>* y, const VectorBase<Scalar>& x, const Scalar& beta )
{
  linear_combination(
    1,Teuchos::arrayArg<Scalar>(Teuchos::ScalarTraits<Scalar>::one())()
    ,Teuchos::arrayArg<const VectorBase<Scalar>*>(&x)()
    ,beta,y
    );
}

template<class Scalar>
void Thyra::V_V( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{
  assign(&*y,x);
}

template<class Scalar>
void Thyra::V_S( VectorBase<Scalar>* y, const Scalar& alpha )
{
  assign(&*y,alpha);
}

template<class Scalar>
void Thyra::V_VpV( VectorBase<Scalar>* z, const VectorBase<Scalar>& x, const VectorBase<Scalar>& y )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination(
    2,Teuchos::arrayArg<Scalar>(ST::one(),ST::one())()
    ,Teuchos::arrayArg<const VectorBase<Scalar>*>(&x,&y)()
    ,ST::zero(),z
    );
}

template<class Scalar>
void Thyra::V_VmV( VectorBase<Scalar>* z, const VectorBase<Scalar>& x, const VectorBase<Scalar>& y )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination(
    2,Teuchos::arrayArg<Scalar>(ST::one(),Scalar(-ST::one()))()
    ,Teuchos::arrayArg<const VectorBase<Scalar>*>(&x,&y)()
    ,ST::zero(),z
    );
}

template<class Scalar>
void Thyra::V_StVpV( VectorBase<Scalar>* z, const Scalar &alpha, const VectorBase<Scalar>& x, const VectorBase<Scalar>& y )
{
  linear_combination(
    2,Teuchos::arrayArg<Scalar>(alpha,Teuchos::ScalarTraits<Scalar>::one())()
    ,Teuchos::arrayArg<const VectorBase<Scalar>*>(&x,&y)()
    ,Teuchos::ScalarTraits<Scalar>::zero(),z
    );
}

template<class Scalar>
void Thyra::V_StVpStV( VectorBase<Scalar>* z, const Scalar &alpha, const VectorBase<Scalar>& x, const Scalar &beta, const VectorBase<Scalar>& y )
{
  linear_combination(
    2,Teuchos::arrayArg<Scalar>(alpha,beta)()
    ,Teuchos::arrayArg<const VectorBase<Scalar>*>(&x,&y)()
    ,Teuchos::ScalarTraits<Scalar>::zero(),z
    );
}

//
// For real types only
//

template<class Scalar>
Scalar Thyra::min( const VectorBase<Scalar>& x ) {
  RTOpPack::ROpMin<Scalar> min_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> min_targ = min_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &x };
  applyOp<Scalar>(min_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*min_targ);
  return min_op(*min_targ);
}

template<class Scalar>
void Thyra::min( const VectorBase<Scalar>& x, Scalar *minEle, Index *minIndex )
{
  RTOpPack::ROpMinIndex<Scalar> min_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> min_targ = min_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &x };
  applyOp<Scalar>(min_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*min_targ);
  RTOpPack::ScalarIndex<Scalar> scalarIndex = min_op(*min_targ);
  *minEle   = scalarIndex.scalar;
  *minIndex = scalarIndex.index;
}

template<class Scalar>
void Thyra::minGreaterThanBound( const VectorBase<Scalar>& x, const Scalar &bound, Scalar *minEle, Index *minIndex )
{
  RTOpPack::ROpMinIndexGreaterThanBound<Scalar> min_op(bound);
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> min_targ = min_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &x };
  applyOp<Scalar>(min_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*min_targ);
  RTOpPack::ScalarIndex<Scalar> scalarIndex = min_op(*min_targ);
  *minEle   = scalarIndex.scalar;
  *minIndex = scalarIndex.index;
}

template<class Scalar>
Scalar Thyra::max( const VectorBase<Scalar>& x )
{
  RTOpPack::ROpMax<Scalar> max_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> max_targ = max_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &x };
  applyOp<Scalar>(max_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*max_targ);
  return max_op(*max_targ);
}

template<class Scalar>
void Thyra::max( const VectorBase<Scalar>& x, Scalar *maxEle, Index *maxIndex )
{
  RTOpPack::ROpMaxIndex<Scalar> max_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> max_targ = max_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &x };
  applyOp<Scalar>(max_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*max_targ);
  RTOpPack::ScalarIndex<Scalar> scalarIndex = max_op(*max_targ);
  *maxEle   = scalarIndex.scalar;
  *maxIndex = scalarIndex.index;
}

template<class Scalar>
void Thyra::maxLessThanBound( const VectorBase<Scalar>& x, const Scalar &bound, Scalar *maxEle, Index *maxIndex )
{
  RTOpPack::ROpMaxIndexLessThanBound<Scalar> max_op(bound);
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> max_targ = max_op.reduct_obj_create();
  const VectorBase<Scalar>* vecs[] = { &x };
  applyOp<Scalar>(max_op,1,vecs,0,static_cast<VectorBase<Scalar>**>(NULL),&*max_targ);
  RTOpPack::ScalarIndex<Scalar> scalarIndex = max_op(*max_targ);
  *maxEle   = scalarIndex.scalar;
  *maxIndex = scalarIndex.index;
}

#endif // THYRA_VECTOR_STD_OPS_HPP



