// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_STD_OPS_DECL_HPP
#define THYRA_MULTI_VECTOR_STD_OPS_DECL_HPP

#include "Thyra_MultiVectorBase.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNorm2.hpp"
#include "RTOpPack_ROpNormInf.hpp"

namespace Thyra {


/** \brief Column-wise multi-vector natural norm.
 *
 * \param V [in]
 *
 * \param norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of the
 * natural norms <tt>dot[j] = sqrt(scalarProd(*V.col(j),*V.col(j)))</tt>, for
 * <tt>j=0...m-1</tt>, computed using a single reduction.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void norms( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms );


/** \brief Column-wise multi-vector reductions.
 *
 * \param V [in]
 *
 * \param normOp [in] A reduction operator consistent with the interface to
 * <tt>RTOpPack::ROpScalarReductionBase</tt> that defines the norm operation.
 *
 * \param norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of
 * one-norms <tt>dot[j] = {some norm}(*V.col(j))</tt>, for <tt>j=0...m-1</tt>,
 * computed using a single reduction.
 *
 * \relates MultiVectorBase
 */
template<class Scalar, class NormOp>
void reductions( const MultiVectorBase<Scalar>& V, const NormOp &op,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms );


/** \brief Column-wise multi-vector one norm.
 *
 * \param V [in]
 *
 * \param norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of
 * one-norms <tt>dot[j] = norm_1(*V.col(j))</tt>, for <tt>j=0...m-1</tt>,
 * computed using a single reduction.
 *
 * This function simply calls <tt>reductions()</tt> using
 * <tt>RTOpPack::ROpNorm1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void norms_1( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms );


/** \brief Column-wise multi-vector 2 (Euclidean) norm.
 *
 * \param V [in]
 *
 * \param norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of
 * one-norms <tt>dot[j] = norm_2(*V.col(j))</tt>, for <tt>j=0...m-1</tt>,
 * computed using a single reduction.
 *
 * This function simply calls <tt>reductions()</tt> using
 * <tt>RTOpPack::ROpNorm2</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void norms_2( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms );


/** \brief Column-wise multi-vector infinity norm.
 *
 * \param V [in]
 *
 * \param norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of
 * one-norms <tt>dot[j] = norm_inf(*V.col(j))</tt>, for <tt>j=0...m-1</tt>,
 * computed using a single reduction.
 *
 * This function simply calls <tt>reductions()</tt> using
 * <tt>RTOpPack::ROpNormInf</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void norms_inf( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms );


/** \brief Column-wise multi-vector infinity norm.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
Array<typename ScalarTraits<Scalar>::magnitudeType>
norms_inf( const MultiVectorBase<Scalar>& V );


/** \brief Multi-vector dot product.
 *
 * \param V1 [in]
 *
 * \param V2 [in]
 *
 * \param dots [out] Array (size <tt>m = V1->domain()->dim()</tt>) of the dot
 * products <tt>dot[j] = dot(*V1.col(j),*V2.col(j))</tt>, for
 * <tt>j=0...m-1</tt>, computed using a single reduction.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void dots( const MultiVectorBase<Scalar>& V1, const MultiVectorBase<Scalar>& V2,
  const ArrayView<Scalar> &dots );


/** \brief Multi-vector column sum
 *
 * \param V [in]
 *
 * \param sums [outt] Array (size <tt>m = V->domain()->dim()</tt>) of the sums
 * products <tt>sum[j] = sum(*V.col(j))</tt>, for <tt>j=0...m-1</tt>, computed
 * using a single reduction.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void sums( const MultiVectorBase<Scalar>& V, const ArrayView<Scalar> &sums );


/** \brief Take the induced matrix one norm of a multi-vector.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
typename ScalarTraits<Scalar>::magnitudeType
norm_1( const MultiVectorBase<Scalar>& V );


/** \brief V = alpha*V.
 *
 * Note, if alpha==0.0 then V=0.0 is performed, and if alpha==1.0 then nothing
 * is done.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void scale( Scalar alpha, const Ptr<MultiVectorBase<Scalar> > &V );


/** \brief A*U + V -> V (where A is a diagonal matrix with diagonal a).
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void scaleUpdate( const VectorBase<Scalar>& a, const MultiVectorBase<Scalar>& U,
  const Ptr<MultiVectorBase<Scalar> > &V );

/** \brief V = alpha.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void assign( const Ptr<MultiVectorBase<Scalar> > &V, Scalar alpha );

/** \brief V = U.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void assign( const Ptr<MultiVectorBase<Scalar> > &V,
  const MultiVectorBase<Scalar>& U );


/** \brief alpha*U + V -> V.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void update( Scalar alpha, const MultiVectorBase<Scalar>& U,
  const Ptr<MultiVectorBase<Scalar> > &V );


/** \brief alpha[j]*beta*U(j) + V(j) - > V(j), for j = 0 ,,,
 *
 * \relates MultiVectorBase
 * U.domain()->dim()-1.
 */
template<class Scalar>
void update(
  const ArrayView<const Scalar> &alpha,
  Scalar beta,
  const MultiVectorBase<Scalar>& U,
  const Ptr<MultiVectorBase<Scalar> > &V
  );


/** \brief U(j) + alpha[j]*beta*V(j) - > V(j), for j = 0 ,,,
 * U.domain()->dim()-1.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void update(
  const MultiVectorBase<Scalar>& U,
  const ArrayView<const Scalar> &alpha,
  Scalar beta,
  const Ptr<MultiVectorBase<Scalar> > &V
  );


/** \brief <tt>Y.col(j)(i) = beta*Y.col(j)(i) + sum( alpha[k]*X[k].col(j)(i),
 * k=0...m-1 )</tt>, for <tt>i = 0...Y->range()->dim()-1</tt>, <tt>j =
 * 0...Y->domain()->dim()-1</tt>.
 *
 * \param alpha [in] Array (length <tt>m</tt>) of input scalars.
 *
 * \param X [in] Array (length <tt>m</tt>) of input multi-vectors.
 *
 * \param beta [in] Scalar multiplier for Y
 *
 * \param Y [in/out] Target multi-vector that is the result of the linear
 * combination.
 *
 * This function implements a general linear combination:
 \verbatim
 Y.col(j)(i) = beta*Y.col(j)(i) + alpha[0]*X[0].col(j)(i) + alpha[1]*X[1].col(j)(i) + ... + alpha[m-1]*X[m-1].col(j)(i)

    for:
        i = 0...y->space()->dim()-1
        j = 0...y->domain()->dim()-1

 \endverbatim
 * and does so on a single call to <tt>MultiVectorBase::applyOp()</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void linear_combination(
  const ArrayView<const Scalar> &alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &X,
  const Scalar &beta,
  const Ptr<MultiVectorBase<Scalar> > &Y
  );


/** \brief Generate a random multi-vector with elements uniformly distributed
 * elements.
 * 
 * The elements <tt>get_ele(*V->col(j))</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using <tt>seed_randomize()</tt>
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, const Ptr<MultiVectorBase<Scalar> > &V );


/** \brief <tt>Z(i,j) *= alpha, i = 0...Z->range()->dim()-1, j =
 * 0...Z->domain()->dim()-1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void Vt_S( const Ptr<MultiVectorBase<Scalar> > &Z, const Scalar& alpha );


/** \brief <tt>Z(i,j) += alpha, i = 0...Z->range()->dim()-1, j =
 * 0...Z->domain()->dim()-1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void Vp_S( const Ptr<MultiVectorBase<Scalar> > &Z, const Scalar& alpha );


/** \brief <tt>Z(i,j) += X(i,j), i = 0...Z->range()->dim()-1, j =
 * 0...Z->domain()->dim()-1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void Vp_V( const Ptr<MultiVectorBase<Scalar> > &Z,
  const MultiVectorBase<Scalar>& X );


/** \brief <tt>Z(i,j) = X(i,j) + Y(i,j), i = 0...Z->range()->dim()-1, j =
 * 0...Z->domain()->dim()-1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void V_VpV( const Ptr<MultiVectorBase<Scalar> > &Z,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y );


/** \brief <tt>Z(i,j) = X(i,j) - Y(i,j), i = 0...Z->range()->dim()-1, j =
 * 0...Z->domain()->dim()-1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void V_VmV( const Ptr<MultiVectorBase<Scalar> > &Z,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y );


/** \brief <tt>Z(i,j) = alpha*X(i,j) + Y(i), i = 0...z->space()->dim()-1</tt>,
 * , j = 0...Z->domain()->dim()-1</tt>.
 *
 * \relates MultiVectorBase
 */
template<class Scalar>
void V_StVpV( const Ptr<MultiVectorBase<Scalar> > &Z, const Scalar &alpha,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y );

#ifndef THYRA_HIDE_DEPRECATED_CODE
/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void norms( const MultiVectorBase<Scalar>& V,
  typename ScalarTraits<Scalar>::magnitudeType norms_out[] )
{ norms(V, Teuchos::arrayView(norms_out, V.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar, class NormOp>
THYRA_DEPRECATED
void reductions( const MultiVectorBase<Scalar>& V, const NormOp &op,
  typename ScalarTraits<Scalar>::magnitudeType norms_out[] )
{ reductions(V, op, Teuchos::arrayView(norms_out, V.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void norms_1( const MultiVectorBase<Scalar>& V,
  typename ScalarTraits<Scalar>::magnitudeType norms_out[] )
{ norms_1(V, Teuchos::arrayView(norms_out, V.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void norms_2( const MultiVectorBase<Scalar>& V,
  typename ScalarTraits<Scalar>::magnitudeType norms_out[] )
{ norms_2(V, Teuchos::arrayView(norms_out, V.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void norms_inf( const MultiVectorBase<Scalar>& V,
  typename ScalarTraits<Scalar>::magnitudeType norms_out[] )
{ norms_inf<Scalar>(V, Teuchos::arrayView(norms_out, V.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void dots( const MultiVectorBase<Scalar>& V1, const MultiVectorBase<Scalar>& V2,
  Scalar dots_out[] )
{ dots<Scalar>(V1, V2, Teuchos::arrayView(dots_out, V1.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void sums( const MultiVectorBase<Scalar>& V, Scalar sums_out[] )
{ sums<Scalar>(V, Teuchos::arrayView(sums_out, V.domain()->dim())); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void scale( Scalar alpha, MultiVectorBase<Scalar>* V )
{ scale(alpha, Teuchos::ptr(V)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void scaleUpdate( const VectorBase<Scalar>& a, const MultiVectorBase<Scalar>& U,
  MultiVectorBase<Scalar>* V )
{ scaleUpdate(a, U, Teuchos::ptr(V)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void assign( MultiVectorBase<Scalar>* V, Scalar alpha )
{ assign(Teuchos::ptr(V), alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void assign( MultiVectorBase<Scalar>* V, const MultiVectorBase<Scalar>& U )
{ assign(Teuchos::ptr(V), U); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void update( Scalar alpha, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V )
{ update(alpha, U, Teuchos::ptr(V)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void update( const Scalar alpha[], Scalar beta, const MultiVectorBase<Scalar>& U,
  MultiVectorBase<Scalar>* V )
{ update(Teuchos::arrayView(alpha, U.domain()->dim()), beta, U, Teuchos::ptr(V)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void update( const MultiVectorBase<Scalar>& U, Scalar alpha[], Scalar beta,
  MultiVectorBase<Scalar>* V )
{ update(U, Teuchos::arrayView(alpha, U.domain()->dim()), beta, Teuchos::ptr(V)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void linear_combination(
  const int m
  ,const Scalar alpha[]
  ,const MultiVectorBase<Scalar>* X_in[]
  ,const Scalar &beta
  ,MultiVectorBase<Scalar> *Y
  )
{
  Array<Ptr<const MultiVectorBase<Scalar> > > X(m);
  for ( int k = 0; k < m; ++k )
    X[k] = Teuchos::ptr(X_in[k]);
  linear_combination<Scalar>(
    Teuchos::arrayView(alpha,m), X(), beta, Teuchos::ptr(Y) );
}


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void randomize( Scalar l, Scalar u, MultiVectorBase<Scalar>* V )
{ randomize(l, u, Teuchos::ptr(V)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vt_S( MultiVectorBase<Scalar>* Z, const Scalar& alpha )
{ Vt_S(Teuchos::ptr(Z), alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vp_S( MultiVectorBase<Scalar>* Z, const Scalar& alpha )
{ Vp_S(Teuchos::ptr(Z), alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vp_V( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X )
{ Vp_V(Teuchos::ptr(Z), X); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_VpV( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X,
  const MultiVectorBase<Scalar>& Y )
{ V_VpV(Teuchos::ptr(Z), X, Y); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_VmV( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y )
{ V_VmV(Teuchos::ptr(Z), X, Y); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_StVpV(
  MultiVectorBase<Scalar>* Z, const Scalar &alpha,
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y 
  )
{ V_StVpV(Teuchos::ptr(Z), alpha, X, Y); }
#endif // THYRA_HIDE_DEPRECATED_CODE

} // end namespace Thyra


// /////////////////////////////////////
// Inline functions


template<class Scalar>
inline
void Thyra::norms_1( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms )
{
  reductions<Scalar>(V, RTOpPack::ROpNorm1<Scalar>(), norms);
}


template<class Scalar>
inline
void Thyra::norms_2( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms )
{
  reductions<Scalar>(V, RTOpPack::ROpNorm2<Scalar>(), norms);
}


template<class Scalar>
inline
void Thyra::norms_inf( const MultiVectorBase<Scalar>& V,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms )
{
  reductions<Scalar>(V, RTOpPack::ROpNormInf<Scalar>(), norms);
}


template<class Scalar>
Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
Thyra::norms_inf( const MultiVectorBase<Scalar>& V )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  Array<ScalarMag> norms(V.domain()->dim());
  Thyra::norms_inf<Scalar>(V, norms());
  return norms;
}


// /////////////////////////////////////////////
// Other implementations


template<class Scalar, class NormOp>
void Thyra::reductions( const MultiVectorBase<Scalar>& V, const NormOp &op,
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  const int m = V.domain()->dim();
  Array<RCP<RTOpPack::ReductTarget> > rcp_op_targs(m);
  Array<Ptr<RTOpPack::ReductTarget> > op_targs(m);
  for( int kc = 0; kc < m; ++kc ) {
    rcp_op_targs[kc] = op.reduct_obj_create();
    op_targs[kc] = rcp_op_targs[kc].ptr();
  }
  applyOp<Scalar>(op, tuple(ptrInArg(V)),
    ArrayView<Ptr<MultiVectorBase<Scalar> > >(null),
    op_targs );
  for( int kc = 0; kc < m; ++kc ) {
    norms[kc] = op(*op_targs[kc]);
  }
}


#endif // THYRA_MULTI_VECTOR_STD_OPS_DECL_HPP
