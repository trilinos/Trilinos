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


#ifndef THYRA_VECTOR_STD_OPS_DECL_HPP
#define THYRA_VECTOR_STD_OPS_DECL_HPP


#include "Thyra_OperatorVectorTypes.hpp"


namespace Thyra {


/** \brief Sum of vector elements:
 * <tt>result = sum( v(i), i = 0...v.space()->dim()-1 )</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar sum( const VectorBase<Scalar>& v );


/** \brief Scalar product <tt>result = <x,y></tt>.
 *
 * Returns <tt>x.space()->scalarProd(x,y)</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y );


/** \brief Inner/Scalar product <tt>result = <x,y></tt>.
 *
 * Returns <tt>x.space()->scalarProd(x,y)</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar inner( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y );


/** \brief Natural norm: <tt>result = sqrt(<v,v>)</tt>.
 *
 * Returns
 * <tt>Teuchos::ScalarTraits<Scalar>::squareroot(v.space()->scalarProd(v,v))</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm( const VectorBase<Scalar>& v );


/** \brief One (1) norm: <tt>result = ||v||1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_1( const VectorBase<Scalar>& v );


/** \brief Euclidean (2) norm: <tt>result = ||v||2</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_2( const VectorBase<Scalar>& v );


/** \brief Weighted Euclidean (2) norm:
 * <tt>result = sqrt( sum( w(i)*conj(v(i))*v(i)) )</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_2( const VectorBase<Scalar> &w, const VectorBase<Scalar>& v );


/** \brief Infinity norm: <tt>result = ||v||inf</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_inf( const VectorBase<Scalar>& v_rhs );


/** \brief Dot product: <tt>result = conj(x)'*y</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar dot( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y );


/** \brief Get single element: <tt>result = v(i)</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar get_ele( const VectorBase<Scalar>& v, Ordinal i );


/** \brief Set single element: <tt>v(i) = alpha</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void set_ele( Ordinal i, Scalar alpha, const Ptr<VectorBase<Scalar> > &v );


/** \brief Assign all elements to a scalar:
 * <tt>y(i) = alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void put_scalar( const Scalar& alpha, const Ptr<VectorBase<Scalar> > &y );


/** \brief Vector assignment:
 * <tt>y(i) = x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void copy( const VectorBase<Scalar>& x, const Ptr<VectorBase<Scalar> > &y );


/** \brief Add a scalar to all elements:
 * <tt>y(i) += alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void add_scalar( const Scalar& alpha, const Ptr<VectorBase<Scalar> > &y );


/** \brief Scale all elements by a scalar:
 * <tt>y(i) *= alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>y = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 *
 * \relates VectorBase
 */
template<class Scalar>
void scale( const Scalar& alpha, const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise absolute value:
 * <tt>y(i) = abs(x(i)), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void abs( const VectorBase<Scalar> &x, const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise reciprocal:
 * <tt>y(i) = 1/x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void reciprocal( const VectorBase<Scalar> &x, const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise product update:
 * <tt>y(i) += alpha * x(i) * v(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void ele_wise_prod( const Scalar& alpha, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& v, const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise conjugate product update:
 * <tt>y(i) += alpha * conj(x(i)) * v(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void ele_wise_conj_prod( const Scalar& alpha, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& v, const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise scaling:
 * <tt>y(i) *= x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void ele_wise_scale( const VectorBase<Scalar>& x, const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise product update:
 * <tt>y(i) += alpha * x(i) * v(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void Vp_StVtV( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha,
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& v);


/** \brief Element-wise product update:
 * <tt>y(i) *= alpha * x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void ele_wise_prod_update( const Scalar& alpha, const VectorBase<Scalar>& x,
  const Ptr<VectorBase<Scalar> > &y );


/** \brief Element-wise product update:
 * <tt>y(i) *= alpha * x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void Vt_StV( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha,
  const VectorBase<Scalar> &x );


/** \brief Element-wise division update:
 * <tt>y(i) += alpha * x(i) / v(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void ele_wise_divide( const Scalar& alpha, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& v, const Ptr<VectorBase<Scalar> > &y );


/** \brief Linear combination:
 * <tt>y(i) = beta*y(i) + sum( alpha[k]*x[k](i), k=0...m-1 ), i = 0...y->space()->dim()-1</tt>.
 *
 * \param m [in] Number of vectors x[]
 *
 * \param alpha [in] Array (length <tt>m</tt>) of input scalars.
 *
 * \param x [in] Array (length <tt>m</tt>) of input vectors.
 *
 * \param beta [in] Scalar multiplier for y
 *
 * \param y [in/out] Target vector that is the result of the linear
 * combination.
 *
 * This function implements a general linear combination:
 \verbatim
   y(i) = beta*y(i) + alpha[0]*x[0](i) + alpha[1]*x[1](i)
          + ... + alpha[m-1]*x[m-1](i), i = 0...y->space()->dim()-1
 \endverbatim
 *
 * \relates VectorBase
 */
template<class Scalar>
void linear_combination(
  const ArrayView<const Scalar> &alpha,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &x,
  const Scalar &beta,
  const Ptr<VectorBase<Scalar> > &y
  );


/** \brief Seed the random number generator used in <tt>randomize()</tt>.
 *
 * \param s [in] The seed for the random number generator.
 *
 * Note, this just calls
 * <tt>Teuchos::TOpRandomize<Scalar>::set_static_seed(s)</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void seed_randomize( unsigned int s );


/** \brief Random vector generation:
 * <tt>v(i) = rand(l,u), , i = 1...v->space()->dim()</tt>.
 * 
 * The elements <tt>v->getEle(i)</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using the above <tt>seed_randomize()</tt> function.
 *
 * \relates VectorBase
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, const Ptr<VectorBase<Scalar> > &v );


/** \brief Assign all elements to a scalar:
 * <tt>y(i) = alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void assign( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha );


/** \brief Vector assignment:
 * <tt>y(i) = x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void assign( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x );


/** \brief Add a scalar to all elements:
 * <tt>y(i) += alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void Vp_S( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha );


/** \brief Scale all elements by a scalar:
 * <tt>y(i) *= alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>y = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 *
 * \relates VectorBase
 */
template<class Scalar>
void Vt_S( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha );


/** \brief Assign scaled vector:
 * <tt>y(i) = alpha * x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_StV( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha,
  const VectorBase<Scalar> &x );


/** \brief AXPY:
 * <tt>y(i) = alpha * x(i) + y(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void Vp_StV( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha,
  const VectorBase<Scalar>& x );


/** \brief <tt>y(i) = x(i) + beta*y(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void Vp_V(
  const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x,
  const Scalar& beta = static_cast<Scalar>(1.0)
  );


/** \brief <tt>y(i) = x(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_V( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar>& x );


/** \brief <tt>y(i) = alpha, i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_S( const Ptr<VectorBase<Scalar> > &y, const Scalar& alpha );


/** \brief <tt>z(i) = x(i) + y(i), i = 0...z->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_VpV( const Ptr<VectorBase<Scalar> > &z, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& y );


/** \brief <tt>z(i) = x(i) - y(i), i = 0...z->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_VmV( const Ptr<VectorBase<Scalar> > &z, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& y );


/** \brief <tt>z(i) = alpha*x(i) + y(i), i = 0...z->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_StVpV( const Ptr<VectorBase<Scalar> > &z, const Scalar &alpha,
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& y );


/** \brief <tt>z(i) = x(i) + alpha*y(i), i = 0...z->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_VpStV( const Ptr<VectorBase<Scalar> > &z,
  const VectorBase<Scalar>& x,
  const Scalar &alpha, const VectorBase<Scalar>& y );


/** \brief <tt>z(i) = alpha*x(i) + beta*y(i), i = 0...z->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void V_StVpStV( const Ptr<VectorBase<Scalar> > &z, const Scalar &alpha,
  const VectorBase<Scalar>& x, const Scalar &beta, const VectorBase<Scalar>& y );


/** \brief Min element: <tt>result = min{ x(i), i = 0...x.space()->dim()-1 } </tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar min( const VectorBase<Scalar>& x );


/** \brief Min element and its index: Returns <tt>maxEle = x(k)</tt>
 * and <tt>maxIndex = k</tt> such that <tt>x(k) <= x(i)</tt> for all
 * <tt>i = 0...x.space()->dim()-1</tt>.
 *
 * \param x [in] Input vector.
 *
 * \param minEle [out] The minimum element value.
 *
 * \param maxindex [out] The global index of the minimum element.  If there is
 * more than one element with the maximum entry then this returns the lowest
 * index in order to make the output independent of the order of operations.
 *
 * Preconditions:<ul>
 * <li><tt>minEle!=NULL</tt>
 * <li><tt>minIndex!=NULL</tt>
 * </ul>
 *
 * \relates VectorBase
 */
template<class Scalar>
void min( const VectorBase<Scalar>& x,
  const Ptr<Scalar> &maxEle, const Ptr<Ordinal> &maxIndex );


/** \brief Minimum element greater than some bound and its index:
 * Returns <tt>minEle = x(k)</tt> and <tt>minIndex = k</tt> such that
 * <tt>x(k) <= x(i)</tt> for all <tt>i</tt> where <tt>x(i) >
 * bound</tt>.
 *
 * \param x [in] Input vector.
 *
 * \param bound [in] The upper bound
 *
 * \param minEle [out] The minimum element value as defined above.
 *
 * \param minIndex [out] The global index of the maximum element.  If there is
 * more than one element with the minimum value then this returns the lowest
 * index in order to make the output independent of the order of operations.
 * If no entries are less than <tt>bound</tt> then <tt>minIndex < 0</tt> on
 * return.
 *
 * Preconditions:<ul>
 * <li><tt>minEle!=NULL</tt>
 * <li><tt>minIndex!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li>If <tt>*minIndex > 0</tt> then such an element was found.
 * <li>If <tt>*minIndex < 0</tt> then no such element was found.
 * </ul>
 *
 * \relates VectorBase
 */
template<class Scalar>
void minGreaterThanBound( const VectorBase<Scalar>& x, const Scalar &bound,
  const Ptr<Scalar> &minEle, const Ptr<Ordinal> &minIndex );


/** \brief Max element: <tt>result = max{ x(i), i = 1...n } </tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
Scalar max( const VectorBase<Scalar>& x );


/** \brief Max element and its index: Returns <tt>maxEle = x(k)</tt>
 * and <tt>maxIndex = k</tt> such that <tt>x(k) >= x(i)</tt> for
 * <tt>i = 0...x.space()->dim()-1</tt>.
 *
 * \param x [in] Input vector.
 *
 * \param maxEle [out] The maximum element value.
 *
 * \param maxindex [out] The global index of the maximum element.  If there is
 * more than one element with the maximum value then this returns the lowest
 * index in order to make the output independent of the order of operations.
 *
 * Preconditions:<ul>
 * <li><tt>maxEle!=NULL</tt>
 * <li><tt>maxIndex!=NULL</tt>
 * </ul>
 *
 * \relates VectorBase
 */
template<class Scalar>
void max( const VectorBase<Scalar>& x,
  const Ptr<Scalar> &maxEle, const Ptr<Ordinal> &maxIndex );


/** \brief Max element less than bound and its index: Returns <tt>maxEle =
 * x(k)</tt> and <tt>maxIndex = k</tt> such that <tt>x(k) >= x(i)</tt> for all
 * <tt>i</tt> where <tt>x(i) < bound</tt>.
 *
 * \param x [in] Input vector.
 *
 * \param bound [in] The upper bound
 *
 * \param maxEle [out] The maximum element value as defined above.
 *
 * \param maxindex [out] The global index of the maximum element.  If there is
 * more than one element with the maximum index then this returns the lowest
 * index in order to make the output independent of the order of operations.
 * If no entries are less than <tt>bound</tt> then <tt>minIndex < 0</tt> on
 * return.
 *
 * Preconditions:<ul>
 * <li><tt>maxEle!=NULL</tt>
 * <li><tt>maxIndex!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li>If <tt>*maxIndex > 0</tt> then such an element was found.
 * <li>If <tt>*maxIndex < 0</tt> then no such element was found.
 * </ul>
 *
 * \relates VectorBase
 */
template<class Scalar>
void maxLessThanBound( const VectorBase<Scalar>& x, const Scalar &bound,
  const Ptr<Scalar> &maxEle, const Ptr<Ordinal> &maxIndex );


#ifndef THYRA_HIDE_DEPRECATED_CODE


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void set_ele( Ordinal i, Scalar alpha, VectorBase<Scalar>* v )
{ set_ele(i, alpha, Teuchos::ptr(v)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void put_scalar( const Scalar& alpha, VectorBase<Scalar>* y )
{ put_scalar<Scalar>(alpha,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void copy( const VectorBase<Scalar>& x, VectorBase<Scalar>* y )
{ copy(x,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void add_scalar( const Scalar& alpha, VectorBase<Scalar>* y )
{ add_scalar(alpha,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void scale( const Scalar& alpha, VectorBase<Scalar>* y )
{ scale(alpha,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void abs( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{ abs(Teuchos::ptr(y),x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void abs( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar> &x )
{
  abs<Scalar>( x, y );
}


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void reciprocal( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{ reciprocal(Teuchos::ptr(y),x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void reciprocal( const Ptr<VectorBase<Scalar> > &y, const VectorBase<Scalar> &x )
{
  reciprocal<Scalar>( x, y );
}


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void ele_wise_prod( const Scalar& alpha, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& v, VectorBase<Scalar>* y )
{ ele_wise_prod(alpha, x, v, Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void ele_wise_conj_prod( const Scalar& alpha, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& v, VectorBase<Scalar>* y )
{ ele_wise_conj_prod(alpha,x,v,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void Vp_StVtV( VectorBase<Scalar>* y, const Scalar& alpha,
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& v)
{ Vp_StVtV(Teuchos::ptr(y),alpha,x,v); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void ele_wise_prod_update( const Scalar& alpha, const VectorBase<Scalar>& x,
  VectorBase<Scalar>* y )
{ ele_wise_prod_update(alpha,x,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void Vt_StV( VectorBase<Scalar>* y, const Scalar& alpha,
  const VectorBase<Scalar> &x )
{ Vt_StV(Teuchos::ptr(y),alpha,x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void ele_wise_divide( const Scalar& alpha, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& v, VectorBase<Scalar>* y )
{ ele_wise_divide(alpha,x,v,Teuchos::ptr(y)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void linear_combination(
  const int m
  ,const Scalar alpha_in[]
  ,const VectorBase<Scalar>* x_in[]
  ,const Scalar &beta
  ,VectorBase<Scalar> *y
  )
{
  Array<Scalar> alpha(m);
  Array<Ptr<const VectorBase<Scalar> > > x(m);
  for (int k = 0; k < m; ++k) {
    alpha[k] = alpha_in[k];
    x[k] = Teuchos::ptr(x_in[k]);
  }
  linear_combination<Scalar>( alpha, x, beta, Teuchos::ptr(y) );
}


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void randomize( Scalar l, Scalar u, VectorBase<Scalar>* v )
{ randomize(l,u,Teuchos::ptr(v)); }


/** \brief Deprecated. */
template<class Scalar> inline
THYRA_DEPRECATED
void assign( VectorBase<Scalar>* y, const Scalar& alpha )
{ assign(Teuchos::ptr(y),alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void assign( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{ assign(Teuchos::ptr(y), x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vp_S( VectorBase<Scalar>* y, const Scalar& alpha )
{ Vp_S(Teuchos::ptr(y), alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vt_S( VectorBase<Scalar>* y, const Scalar& alpha )
{ Vt_S(Teuchos::ptr(y), alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_StV( VectorBase<Scalar>* y, const Scalar& alpha,
  const VectorBase<Scalar> &x )
{ V_StV(Teuchos::ptr(y), alpha, x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vp_StV( VectorBase<Scalar>* y, const Scalar& alpha,
  const VectorBase<Scalar>& x )
{ Vp_StV(Teuchos::ptr(y), alpha, x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void Vp_V(
  VectorBase<Scalar>* y, const VectorBase<Scalar>& x,
  const Scalar& beta = static_cast<Scalar>(1.0)
  )
{ Vp_V(Teuchos::ptr(y), x, beta); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_V( VectorBase<Scalar>* y, const VectorBase<Scalar>& x )
{ V_V(Teuchos::ptr(y), x); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_S( VectorBase<Scalar>* y, const Scalar& alpha )
{ V_S(Teuchos::ptr(y), alpha); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_VpV( VectorBase<Scalar>* z, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& y )
{ V_VpV(Teuchos::ptr(z), x, y); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_VmV( VectorBase<Scalar>* z, const VectorBase<Scalar>& x,
  const VectorBase<Scalar>& y )
{ V_VmV(Teuchos::ptr(z), x, y); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_StVpV( VectorBase<Scalar>* z, const Scalar &alpha,
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& y )
{ V_StVpV(Teuchos::ptr(z), alpha, x, y); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void V_StVpStV( VectorBase<Scalar>* z, const Scalar &alpha,
  const VectorBase<Scalar>& x, const Scalar &beta, const VectorBase<Scalar>& y )
{ V_StVpStV(Teuchos::ptr(z), alpha, x, beta, y); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void min( const VectorBase<Scalar>& x, Scalar *maxEle, Ordinal *maxIndex )
{ min(x, Teuchos::ptr(maxEle), Teuchos::ptr(maxIndex)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void minGreaterThanBound( const VectorBase<Scalar>& x, const Scalar &bound,
  Scalar *minEle, Ordinal *minIndex )
{ minGreaterThanBound(x, bound, Teuchos::ptr(minEle), Teuchos::ptr(minIndex)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void max( const VectorBase<Scalar>& x, Scalar *maxEle, Ordinal *maxIndex )
{ max(x, Teuchos::ptr(maxEle), Teuchos::ptr(maxIndex)); }


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
void maxLessThanBound( const VectorBase<Scalar>& x, const Scalar &bound,
  Scalar *maxEle, Ordinal *maxIndex )
{ maxLessThanBound(x, bound, Teuchos::ptr(maxEle), Teuchos::ptr(maxIndex)); }


#endif // THYRA_HIDE_DEPRECATED_CODE


} // end namespace Thyra


// /////////////////////////
// Inline functions


template<class Scalar>
inline
Scalar Thyra::scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y )
{
  return x.space()->scalarProd(x, y);
}


template<class Scalar>
inline
Scalar Thyra::inner( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y )
{
  return x.space()->scalarProd(x, y);
}


template<class Scalar>
inline
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm( const VectorBase<Scalar>& v )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ST::magnitude(ST::squareroot(v.space()->scalarProd(v, v)));
}


#endif // THYRA_VECTOR_STD_OPS_DECL_HPP
