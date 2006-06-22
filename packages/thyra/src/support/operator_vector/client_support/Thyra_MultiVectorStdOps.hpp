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

#ifndef THYRA_MULTI_VECTOR_STD_OPS_HPP
#define THYRA_MULTI_VECTOR_STD_OPS_HPP

#include "Thyra_MultiVectorStdOpsDecl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
#include "RTOpPack_TOpAssignVectors.hpp"
#include "RTOpPack_TOpAXPY.hpp"
#include "RTOpPack_TOpLinearCombination.hpp"
#include "RTOpPack_TOpScaleVector.hpp"
#include "Teuchos_TestForException.hpp"

template<class Scalar>
void Thyra::norms( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  const int m = V.domain()->dim();
  std::vector<Scalar> prods(m);
  V.range()->scalarProds(V,V,&prods[0]);
  for( int j = 0; j < m; ++j )
    norms[j] = ST::magnitude(ST::squareroot(prods[j]));
}

template<class Scalar, class NormOp>
void Thyra::reductions( const MultiVectorBase<Scalar>& V, const NormOp &op, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] )
{
  int kc;
  const int m = V.domain()->dim();
  std::vector<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >  rcp_op_targs(m);
  std::vector<RTOpPack::ReductTarget*>                        op_targs(m);
  for( kc = 0; kc < m; ++kc ) {
    rcp_op_targs[kc] = op.reduct_obj_create();
    op_targs[kc] = &*rcp_op_targs[kc];
  }
  const MultiVectorBase<Scalar>* multi_vecs[] = { &V,};
  applyOp<Scalar>(op,1,multi_vecs,0,static_cast<MultiVectorBase<Scalar>**>(NULL),&op_targs[0]);
  for( kc = 0; kc < m; ++kc ) {
    norms[kc] = op(*op_targs[kc]);
  }
}

template<class Scalar>
void Thyra::dots( const MultiVectorBase<Scalar>& V1, const MultiVectorBase<Scalar>& V2, Scalar dots[] )
{
  int kc;
  const int m = V1.domain()->dim();
  RTOpPack::ROpDotProd<Scalar> dot_op;
  std::vector<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >  rcp_dot_targs(m);
  std::vector<RTOpPack::ReductTarget*>                        dot_targs(m);
  for( kc = 0; kc < m; ++kc ) {
    rcp_dot_targs[kc] = dot_op.reduct_obj_create();
    dot_targs[kc] = &*rcp_dot_targs[kc];
  }
  const MultiVectorBase<Scalar>* multi_vecs[] = { &V1, &V2 };
  applyOp<Scalar>(dot_op,2,multi_vecs,0,static_cast<MultiVectorBase<Scalar>**>(NULL),&dot_targs[0]); // Cast required by sun compiler
  for( kc = 0; kc < m; ++kc ) {
    dots[kc] = dot_op(*dot_targs[kc]);
  }
}

template<class Scalar>
void Thyra::sums( const MultiVectorBase<Scalar>& V, Scalar sums[] )
{
  int kc;
  const int m = V.domain()->dim();
  RTOpPack::ROpSum<Scalar> sum_op;
  std::vector<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >  rcp_op_targs(m);
  std::vector<RTOpPack::ReductTarget*>                        op_targs(m);
  for( kc = 0; kc < m; ++kc ) {
    rcp_op_targs[kc] = sum_op.reduct_obj_create();
    op_targs[kc] = &*rcp_op_targs[kc];
  }
  const MultiVectorBase<Scalar>* multi_vecs[] = { &V,};
  applyOp<Scalar>(sum_op,1,multi_vecs,0,static_cast<MultiVectorBase<Scalar>**>(NULL),&op_targs[0]);
  for( kc = 0; kc < m; ++kc ) {
    sums[kc] = sum_op(*op_targs[kc]);
  }
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_1( const MultiVectorBase<Scalar>& V )
{
  // Primary column-wise reduction (sum of absolute values)
  RTOpPack::ROpNorm1<Scalar> sum_abs_op;
  // Secondary reduction (max over all columns = induced norm_1 matrix  norm)
  RTOpPack::ROpNormInf<Scalar> max_op;
  // Reduction object (must be same for both sum_abs and max_targ objects)
  Teuchos::RefCountPtr<RTOpPack::ReductTarget>
    max_targ = max_op.reduct_obj_create();
  // Perform the reductions
  const MultiVectorBase<Scalar>* multi_vecs[] = { &V };
  applyOp<Scalar>(sum_abs_op,max_op,1,multi_vecs,0,static_cast<MultiVectorBase<Scalar>**>(NULL),&*max_targ); // Sun complier requires this cast
  // Return the final value
  return max_op(*max_targ);
}

template<class Scalar>
void Thyra::scale( Scalar alpha, MultiVectorBase<Scalar>* V )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
  if(alpha==ST::zero()) {
    assign( V, ST::zero() );
    return;
  }
  if(alpha==ST::one()) {
    return;
  }
  RTOpPack::TOpScaleVector<Scalar> scale_vector_op(alpha);
  MultiVectorBase<Scalar>* targ_multi_vecs[] = { V };
  applyOp<Scalar>(
    scale_vector_op,0,(const MultiVectorBase<Scalar>**)NULL // The SUN compiler requires these casts!
    ,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL
    );
}

template<class Scalar>
void Thyra::scaleUpdate( const VectorBase<Scalar>& a, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
  bool is_compatible = U.range()->isCompatible(*a.space());
  TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.range()->isCompatible(*a.space())==false");
  is_compatible = U.range()->isCompatible(*V->range());
  TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.range()->isCompatible((V->range())==false ");
  is_compatible = U.domain()->isCompatible(*V->domain());
  TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.domain().isCompatible(V->domain())==false ");
#endif
  const int m = U.domain()->dim();
  for( int j = 0; j < m; ++j ) {
    ele_wise_prod( Scalar(1.0), a, *U.col(j), &*V->col(j) ); 
  }
}

template<class Scalar>
void Thyra::assign( MultiVectorBase<Scalar>* V, Scalar alpha )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
  MultiVectorBase<Scalar>* targ_multi_vecs[] = { V };
  applyOp<Scalar>(
    assign_scalar_op,0,(const MultiVectorBase<Scalar>**)NULL // The SUN compiler requires these casts!
    ,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL
    );
}

template<class Scalar>
void Thyra::assign( MultiVectorBase<Scalar>* V, const MultiVectorBase<Scalar>& U )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAssignVectors<Scalar> assign_vectors_op;
  const MultiVectorBase<Scalar>* multi_vecs[]      = { &U };
  MultiVectorBase<Scalar>*       targ_multi_vecs[] = { V   };
  applyOp<Scalar>(
    assign_vectors_op,1,multi_vecs,1,targ_multi_vecs
    ,(RTOpPack::ReductTarget**)NULL // The SUN compiler requires this cast!
    );
}

template<class Scalar>
void Thyra::update( Scalar alpha, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
#endif
  RTOpPack::TOpAXPY<Scalar> axpy_op(alpha);
  const MultiVectorBase<Scalar>* multi_vecs[]       = { &U };
  MultiVectorBase<Scalar>*       targ_multi_vecs[]  = { V  };
  applyOp<Scalar>(axpy_op,1,multi_vecs,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL); // Sun compiler requires this cast
}

template<class Scalar>
void Thyra::update( Scalar alpha[], Scalar beta, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
  bool is_compatible = U.range()->isCompatible(*V->range());
  TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.range()->isCompatible((V->range())==false ");
  is_compatible = U.domain()->isCompatible(*V->domain());
  TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.domain().isCompatible(V->domain())==false ");
#endif
  const int m = U.domain()->dim();
  for( int j = 0; j < m; ++j )
    Vp_StV( V->col(j).get(), alpha[j]*beta, *U.col(j) );
}


template<class Scalar>
void Thyra::update( const MultiVectorBase<Scalar>& U, Scalar alpha[], Scalar beta, MultiVectorBase<Scalar>* V )
{
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
  bool is_compatible = U.range()->isCompatible(*V->range());
    TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.range()->isCompatible((V->range())==false ");
  is_compatible = U.domain()->isCompatible(*V->domain());
    TEST_FOR_EXCEPTION(
    !is_compatible,Exceptions::IncompatibleVectorSpaces
    ,"update(...), Error, U.domain().isCompatible(V->domain())==false ");
#endif
  const int m = U.domain()->dim();
  for( int j = 0; j < m; ++j ) {
    Vt_S( V->col(j).get(), alpha[j]*beta );
    Vp_StV( V->col(j).get(), 1.0, *U.col(j) );
  }
}

template<class Scalar>
void Thyra::linear_combination(
  const int                         m
  ,const Scalar                     alpha[]
  ,const MultiVectorBase<Scalar>*   X[]
  ,const Scalar                     &beta
  ,MultiVectorBase<Scalar>          *Y
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(Y==NULL,std::logic_error,"linear_combination(...), Error!");
#endif
  if( beta == Teuchos::ScalarTraits<Scalar>::one() && m == 1 ) {
    update( alpha[0], *X[0], Y );
    return;
  }
  else if( m == 0 ) {
    scale( beta, Y );
    return;
  }
  RTOpPack::TOpLinearCombination<Scalar> lin_comb_op(m,alpha,beta);
  MultiVectorBase<Scalar>* targ_multi_vecs[] = { Y };
  Thyra::applyOp<Scalar>(lin_comb_op,m,X,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL);  // Cast returned by sun compiler
}

template<class Scalar>
void Thyra::randomize( Scalar l, Scalar u, MultiVectorBase<Scalar>* V )
{
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"randomize(...), Error!");
#endif
  const int m = V->domain()->dim();
  for( int j = 0; j < m; ++j ) {
    randomize( l, u, V->col(j).get() ); // Todo: call applyOp(...) directly!
  }
}

template<class Scalar>
void Thyra::V_VpV( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination(
    2,Teuchos::arrayArg<Scalar>(ST::one(),ST::one())()
    ,Teuchos::arrayArg<const MultiVectorBase<Scalar>*>(&X,&Y)()
    ,ST::zero(),Z
    );
}

template<class Scalar>
void Thyra::V_VmV( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination(
    2,Teuchos::arrayArg<Scalar>(ST::one(),Scalar(-ST::one()))()
    ,Teuchos::arrayArg<const MultiVectorBase<Scalar>*>(&X,&Y)()
    ,ST::zero(),Z
    );
}

#endif // THYRA_MULTI_VECTOR_STD_OPS_HPP
