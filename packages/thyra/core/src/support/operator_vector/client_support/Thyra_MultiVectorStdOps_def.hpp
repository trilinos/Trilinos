// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_STD_OPS_HPP
#define THYRA_MULTI_VECTOR_STD_OPS_HPP

#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


template<class Scalar>
void Thyra::norms( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  const int m = V.domain()->dim();
  Array<Scalar> prods(m);
  V.range()->scalarProds(V, V, prods());
  for ( int j = 0; j < m; ++j )
    norms[j] = ST::magnitude(ST::squareroot(prods[j]));
}


template<class Scalar>
void Thyra::dots( const MultiVectorBase<Scalar>& V1, const MultiVectorBase<Scalar>& V2,
  const ArrayView<Scalar> &dots )
{
  V2.dots(V1, dots);
}


template<class Scalar>
void Thyra::sums( const MultiVectorBase<Scalar>& V, const ArrayView<Scalar> &sums )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  const int m = V.domain()->dim();
  RTOpPack::ROpSum<Scalar> sum_op;
  Array<RCP<RTOpPack::ReductTarget> > rcp_op_targs(m);
  Array<Ptr<RTOpPack::ReductTarget> > op_targs(m);
  for( int kc = 0; kc < m; ++kc ) {
    rcp_op_targs[kc] = sum_op.reduct_obj_create();
    op_targs[kc] = rcp_op_targs[kc].ptr();
  }
  applyOp<Scalar>(sum_op, tuple(ptrInArg(V)),
    ArrayView<const Ptr<MultiVectorBase<Scalar> > >(null), op_targs);
  for( int kc = 0; kc < m; ++kc ) {
    sums[kc] = sum_op(*op_targs[kc]);
  }
}


template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm_1( const MultiVectorBase<Scalar>& V )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  // Primary column-wise reduction (sum of absolute values)
  RTOpPack::ROpNorm1<Scalar> sum_abs_op;
  // Secondary reduction (max over all columns = induced norm_1 matrix norm)
  RTOpPack::ROpNormInf<Scalar> max_op;
  // Reduction object (must be same for both sum_abs and max_targ objects)
  RCP<RTOpPack::ReductTarget>
    max_targ = max_op.reduct_obj_create();
  // Perform the reductions
  Thyra::applyOp<Scalar>(sum_abs_op, max_op, tuple(ptrInArg(V))(), 
    ArrayView<const Ptr<MultiVectorBase<Scalar> > >(null),
    max_targ.ptr());
  // Return the final value
  return max_op(*max_targ);
}


template<class Scalar>
void Thyra::scale( Scalar alpha, const Ptr<MultiVectorBase<Scalar> > &V )
{
  V->scale(alpha);
}


template<class Scalar>
void Thyra::scaleUpdate( const VectorBase<Scalar>& a,
  const MultiVectorBase<Scalar>& U, const Ptr<MultiVectorBase<Scalar> > &V )
{
#ifdef TEUCHOS_DEBUG
  bool is_compatible = U.range()->isCompatible(*a.space());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, Exceptions::IncompatibleVectorSpaces,
    "update(...), Error, U.range()->isCompatible(*a.space())==false" );
  is_compatible = U.range()->isCompatible(*V->range());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, Exceptions::IncompatibleVectorSpaces,
    "update(...), Error, U.range()->isCompatible((V->range())==false" );
  is_compatible = U.domain()->isCompatible(*V->domain());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, Exceptions::IncompatibleVectorSpaces,
    "update(...), Error, U.domain().isCompatible(V->domain())==false" );
#endif
  const int m = U.domain()->dim();
  for( int j = 0; j < m; ++j ) {
    ele_wise_prod<Scalar>( 1.0, a, *U.col(j), V->col(j).ptr() ); 
  }
}


template<class Scalar>
void Thyra::assign( const Ptr<MultiVectorBase<Scalar> > &V, Scalar alpha )
{
  V->assign(alpha);
}


template<class Scalar>
void Thyra::assign( const Ptr<MultiVectorBase<Scalar> > &V,
  const MultiVectorBase<Scalar>& U )
{
  V->assign(U);
}


template<class Scalar>
void Thyra::update( Scalar alpha, const MultiVectorBase<Scalar>& U,
  const Ptr<MultiVectorBase<Scalar> > &V )
{
  V->update(alpha, U);
}


template<class Scalar>
void Thyra::update( const ArrayView<const Scalar> &alpha, Scalar beta,
  const MultiVectorBase<Scalar>& U, const Ptr<MultiVectorBase<Scalar> > &V )
{
#ifdef TEUCHOS_DEBUG
  bool is_compatible = U.range()->isCompatible(*V->range());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, Exceptions::IncompatibleVectorSpaces,
    "update(...), Error, U.range()->isCompatible((V->range())==false");
  is_compatible = U.domain()->isCompatible(*V->domain());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, Exceptions::IncompatibleVectorSpaces,
    "update(...), Error, U.domain().isCompatible(V->domain())==false");
#endif
  const int m = U.domain()->dim();
  for( int j = 0; j < m; ++j )
    Vp_StV<Scalar>( V->col(j).ptr(), alpha[j]*beta, *U.col(j) );
}


template<class Scalar>
void Thyra::update( const MultiVectorBase<Scalar>& U,
  const ArrayView<const Scalar> &alpha, Scalar beta,
  const Ptr<MultiVectorBase<Scalar> > &V )
{
#ifdef TEUCHOS_DEBUG
  bool is_compatible = U.range()->isCompatible(*V->range());
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_compatible, Exceptions::IncompatibleVectorSpaces,
      "update(...), Error, U.range()->isCompatible((V->range())==false");
    is_compatible = U.domain()->isCompatible(*V->domain());
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_compatible, Exceptions::IncompatibleVectorSpaces,
      "update(...), Error, U.domain().isCompatible(V->domain())==false");
#endif
  const int m = U.domain()->dim();
  for( int j = 0; j < m; ++j ) {
    Vt_S<Scalar>( V->col(j).ptr(), alpha[j]*beta );
    Vp_StV<Scalar>( V->col(j).ptr(), 1.0, *U.col(j) );
  }
}


template<class Scalar>
void Thyra::linear_combination(
  const ArrayView<const Scalar> &alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &X,
  const Scalar &beta,
  const Ptr<MultiVectorBase<Scalar> > &Y
  )
{
  Y->linear_combination(alpha, X, beta);
}


template<class Scalar>
void Thyra::randomize( Scalar l, Scalar u,
  const Ptr<MultiVectorBase<Scalar> > &V )
{
  const int m = V->domain()->dim();
  for( int j = 0; j < m; ++j )
    randomize( l, u, V->col(j).ptr() );
  // Todo: call applyOp(...) directly!
}


template<class Scalar>
void Thyra::Vt_S( const Ptr<MultiVectorBase<Scalar> > &Z,
  const Scalar& alpha )
{
  Z->scale(alpha);
}


template<class Scalar>
void Thyra::Vp_S( const Ptr<MultiVectorBase<Scalar> > &Z,
  const Scalar& alpha )
{
  const int m = Z->domain()->dim();
  for( int j = 0; j < m; ++j )
    Vp_S( Z->col(j).ptr(), alpha );
  // Todo: call applyOp(...) directly!
}


template<class Scalar>
void Thyra::Vp_V( const Ptr<MultiVectorBase<Scalar> > &Z,
  const MultiVectorBase<Scalar>& X )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>( tuple(ST::one()), tuple(ptrInArg(X)),
    ST::one(), Z );
}


template<class Scalar>
void Thyra::V_VpV( const Ptr<MultiVectorBase<Scalar> > &Z,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(ST::one(), ST::one()), tuple(ptrInArg(X), ptrInArg(Y)),
    ST::zero(), Z
    );
}


template<class Scalar>
void Thyra::V_VmV( const Ptr<MultiVectorBase<Scalar> > &Z,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(ST::one(), as<Scalar>(-ST::one())), tuple(ptrInArg(X), ptrInArg(Y)),
    ST::zero(), Z
    );
}


template<class Scalar>
void Thyra::V_StVpV(
  const Ptr<MultiVectorBase<Scalar> > &Z, const Scalar &alpha,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y 
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  linear_combination<Scalar>(
    tuple(alpha, ST::one()),  tuple(ptrInArg(X), ptrInArg(Y)),
    ST::zero(), Z
    );
}


//
// Explicit instant macro
//

#define THYRA_MULTI_VECTOR_STD_OPS_INSTANT(SCALAR) \
   \
  template void norms( const MultiVectorBase<SCALAR >& V, \
    const ArrayView<ScalarTraits<SCALAR >::magnitudeType> &norms ); \
   \
  template void dots( const MultiVectorBase<SCALAR >& V1, const MultiVectorBase<SCALAR >& V2, \
    const ArrayView<SCALAR > &dots ); \
   \
  template void sums( const MultiVectorBase<SCALAR >& V, const ArrayView<SCALAR > &sums ); \
   \
  template Teuchos::ScalarTraits<SCALAR >::magnitudeType \
  norm_1( const MultiVectorBase<SCALAR >& V ); \
   \
  template void scale( SCALAR  alpha, const Ptr<MultiVectorBase<SCALAR > > &V ); \
   \
  template void scaleUpdate( const VectorBase<SCALAR >& a, \
    const MultiVectorBase<SCALAR >& U, const Ptr<MultiVectorBase<SCALAR > > &V ); \
   \
  template void assign( const Ptr<MultiVectorBase<SCALAR > > &V, SCALAR  alpha ); \
   \
  template void assign( const Ptr<MultiVectorBase<SCALAR > > &V, \
    const MultiVectorBase<SCALAR >& U ); \
   \
  template void update( SCALAR  alpha, const MultiVectorBase<SCALAR >& U, \
    const Ptr<MultiVectorBase<SCALAR > > &V ); \
   \
  template void update( const ArrayView<const SCALAR > &alpha, SCALAR  beta, \
    const MultiVectorBase<SCALAR >& U, const Ptr<MultiVectorBase<SCALAR > > &V ); \
   \
  template void update( const MultiVectorBase<SCALAR >& U, \
    const ArrayView<const SCALAR > &alpha, SCALAR  beta, \
    const Ptr<MultiVectorBase<SCALAR > > &V ); \
   \
  template void linear_combination( \
    const ArrayView<const SCALAR > &alpha, \
    const ArrayView<const Ptr<const MultiVectorBase<SCALAR > > > &X, \
    const SCALAR  &beta, \
    const Ptr<MultiVectorBase<SCALAR > > &Y \
    ); \
   \
  template void randomize( SCALAR  l, SCALAR  u, \
    const Ptr<MultiVectorBase<SCALAR > > &V ); \
   \
  template void Vt_S( const Ptr<MultiVectorBase<SCALAR > > &Z, \
    const SCALAR & alpha ); \
   \
  template void Vp_S( const Ptr<MultiVectorBase<SCALAR > > &Z, \
    const SCALAR & alpha ); \
   \
  template void Vp_V( const Ptr<MultiVectorBase<SCALAR > > &Z, \
    const MultiVectorBase<SCALAR >& X ); \
   \
  template void V_VpV( const Ptr<MultiVectorBase<SCALAR > > &Z, \
    const MultiVectorBase<SCALAR >& X, const MultiVectorBase<SCALAR >& Y ); \
   \
  template void V_VmV( const Ptr<MultiVectorBase<SCALAR > > &Z, \
    const MultiVectorBase<SCALAR >& X, const MultiVectorBase<SCALAR >& Y ); \
   \
  template void V_StVpV( \
    const Ptr<MultiVectorBase<SCALAR > > &Z, const SCALAR  &alpha, \
    const MultiVectorBase<SCALAR >& X, const MultiVectorBase<SCALAR >& Y  \
    ); \


#endif // THYRA_MULTI_VECTOR_STD_OPS_HPP
