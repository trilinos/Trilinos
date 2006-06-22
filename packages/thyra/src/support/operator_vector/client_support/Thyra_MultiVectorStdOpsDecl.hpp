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

#ifndef THYRA_MULTI_VECTOR_STD_OPS_DECL_HPP
#define THYRA_MULTI_VECTOR_STD_OPS_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNorm2.hpp"
#include "RTOpPack_ROpNormInf.hpp"

namespace Thyra {


/** \defgroup Thyra_Op_Vec_MultiVectorStdOps_grp Collection of all multi-vector operations
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */

/** \defgroup Thyra_Op_Vec_MultiVectorStdOpsAll_grp Collection of multi-vector operations for all scalar types.
 *
 * \ingroup Thyra_Op_Vec_MultiVectorStdOps_grp
 */

/** \defgroup Thyra_Op_Vec_MultiVectorStdOpsAll_names_grp Collection of standard multi-vector operations with text names.
 *
 * \ingroup Thyra_Op_Vec_MultiVectorStdOpsAll_grp
 */
//@{

/** \brief Column-wise multi-vector natural norm.
 *
 * @param  V     [in]
 * @param  norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of the natural norms
 *               <tt>dot[j] = sqrt(scalarProd(*V.col(j),*V.col(j)))</tt>, for <tt>j=0...m-1</tt>,
 *               computed using a single reduction.
 */
template<class Scalar>
void norms( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] );

/** \brief Column-wise multi-vector reductions.
 *
 * @param  V       [in]
 * @param  normOp  [in] A reduction operator consistent with the interface to
 *                 <tt>RTOpPack::ROpScalarReductionBase</tt> that defines the
 *                 norm operation.
 * @param  norms   [out] Array (size <tt>m = V1->domain()->dim()</tt>) of one-norms
 *                 <tt>dot[j] = {some norm}(*V.col(j))</tt>, for <tt>j=0...m-1</tt>, computed using
 *                 a single reduction.
 */
template<class Scalar, class NormOp>
void reductions( const MultiVectorBase<Scalar>& V, const NormOp &op, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] );

/** \brief Column-wise multi-vector one norm.
 *
 * @param  V     [in]
 * @param  norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of one-norms
 *               <tt>dot[j] = norm_1(*V.col(j))</tt>, for <tt>j=0...m-1</tt>, computed using
 *               a single reduction.
 *
 * This function simply calls <tt>reductions()</tt> using <tt>RTOpPack::ROpNorm1</tt>.
 */
template<class Scalar>
void norms_1( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] );

/** \brief Column-wise multi-vector 2 (Euclidean) norm.
 *
 * @param  V     [in]
 * @param  norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of one-norms
 *               <tt>dot[j] = norm_2(*V.col(j))</tt>, for <tt>j=0...m-1</tt>, computed using
 *               a single reduction.
 *
 * This function simply calls <tt>reductions()</tt> using <tt>RTOpPack::ROpNorm2</tt>.
 */
template<class Scalar>
void norms_2( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] );

/** \brief Column-wise multi-vector infinity norm.
 *
 * @param  V     [in]
 * @param  norms [out] Array (size <tt>m = V1->domain()->dim()</tt>) of one-norms
 *               <tt>dot[j] = norm_inf(*V.col(j))</tt>, for <tt>j=0...m-1</tt>,
 *               computed using a single reduction.
 *
 * This function simply calls <tt>reductions()</tt> using <tt>RTOpPack::ROpNormInf</tt>.
 */
template<class Scalar>
void norms_inf( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] );

/** \brief Multi-vector dot product.
 *
 * @param  V1   [in]
 * @param  V2   [in]
 * @param  dots [out] Array (size <tt>m = V1->domain()->dim()</tt>) of the dot products
 *              <tt>dot[j] = dot(*V1.col(j),*V2.col(j))</tt>, for <tt>j=0...m-1</tt>,
 *              computed using a single reduction.
 */
template<class Scalar>
void dots( const MultiVectorBase<Scalar>& V1, const MultiVectorBase<Scalar>& V2, Scalar dots[] );

/** \brief Multi-vector column sum
 *
 * @param  V    [in]
 * @param  sums [outt] Array (size <tt>m = V->domain()->dim()</tt>) of the sums products
 *              <tt>sum[j] = sum(*V.col(j))</tt>, for <tt>j=0...m-1</tt>,
 *              computed using a single reduction.
 */
template<class Scalar>
void sums( const MultiVectorBase<Scalar>& V, Scalar sums[] );

/** \brief Take the induced matrix one norm of a multi-vector.
 *
 * @param V  [in] Input multi-vector
 *
 *Returns a scalar.
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm_1( const MultiVectorBase<Scalar>& V );

/** \brief V = alpha*V.
 *
 * Note, if alpha==0.0, then V=alpha is performed and if alpha==1.0,
 * then nothing is done.
 */
template<class Scalar>
void scale( Scalar alpha, MultiVectorBase<Scalar>* V );

/** \brief A*U + V -> V (where A is a diagonal matrix with diagonal a).
 */
template<class Scalar>
void scaleUpdate( const VectorBase<Scalar>& a, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V );

/** \brief V = alpha.
 */
template<class Scalar>
void assign( MultiVectorBase<Scalar>* V, Scalar alpha );

/** \brief V = U
 */
template<class Scalar>
void assign( MultiVectorBase<Scalar>* V, const MultiVectorBase<Scalar>& U );

/** \brief alpha*U + V -> V
 */
template<class Scalar>
void update( Scalar alpha, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V );

/** \brief alpha[j-1]*beta*U(j) + V(j) - > V(j), for j = 0 ... U.domain()->dim()-1
 */
template<class Scalar>
void update( Scalar alpha[], Scalar beta, const MultiVectorBase<Scalar>& U, MultiVectorBase<Scalar>* V );

/** \brief U(j) + alpha[j-1]*beta*V(j) - > V(j), for j = 0 ... U.domain()->dim()-1
 */
template<class Scalar>
void update( const MultiVectorBase<Scalar>& U, Scalar alpha[], Scalar beta, MultiVectorBase<Scalar>* V );

/** \brief <tt>Y.col(j)(i) = beta*Y.col(j)(i) + sum( alpha[k]*X[k].col(j)(i), k=0...m-1 )</tt>,
 * for <tt>i = 0...Y->range()->dim()-1</tt>, <tt>j = 0...Y->domain()->dim()-1</tt>.
 *
 * @param  m          [in] Number of multi-vectors in X[]
 * @param  alpha      [in] Array (length <tt>m</tt>) of input scalars.
 * @param  X          [in] Array (length <tt>m</tt>) of input multi-vectors.
 * @param  beta       [in] Scalar multiplier for Y
 * @param  Y          [in/out] Target multi-vector that is the result of the linear combination.
 *
 * This function implements a general linear combination:
 \verbatim
 Y.col(j)(i) = beta*Y.col(j)(i) + alpha[0]*X[0].col(j)(i) + alpha[1]*X[1].col(j)(i) + ... + alpha[m-1]*X[m-1].col(j)(i)

    for:
        i = 0...y->space()->dim()-1
        j = 0...y->domain()->dim()-1

 \endverbatim
 * and does so on a single call to <tt>MultiVectorBase::applyOp()</tt>.
 */
template<class Scalar>
void linear_combination(
  const int                         m
  ,const Scalar                     alpha[]
  ,const MultiVectorBase<Scalar>*   X[]
  ,const Scalar                     &beta
  ,MultiVectorBase<Scalar>          *Y
  );

/** \brief Generate a random multi-vector with elements uniformly distributed
 * elements.
 * 
 * The elements <tt>get_ele(*V->col(j))</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using <tt>seed_randomize()</tt>
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, MultiVectorBase<Scalar>* V );

//@}

/** \defgroup Thyra_Op_Vec_MultiVectorStdOpsAll_LA_names_grp Collection of standard multi-vector operations using linear algebra naming convention.
 *
 * These functions a just simpler ways to call the functions defined
 * \ref Thyra_Op_Vec_MultiVectorStdOpsAll_names_grp "here".
 *
 * The convention used here is described in the short note <a
 * href="./LinearAlgebraFunctionConvention.pdf">A Simple Convention for the
 * Specification of Linear Algebra Function Prototypes in C++ </a>.
 *
 * \ingroup Thyra_Op_Vec_MultiVectorStdOpsAll_grp
 */
//@{

/** \brief <tt>Z(i,j) = X(i,j) + Y(i,j), i = 0...Z->range()->dim()-1, j = 0...Z->domain()->dim()-1</tt>.
 */
template<class Scalar>
void V_VpV( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y );

/** \brief <tt>Z(i,j) = X(i,j) - Y(i,j), i = 0...Z->range()->dim()-1, j = 0...Z->domain()->dim()-1</tt>.
 */
template<class Scalar>
void V_VmV( MultiVectorBase<Scalar>* Z, const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y );

//@}

} // end namespace Thyra

// /////////////////////////////////////
// Inline functions

template<class Scalar>
inline
void Thyra::norms_1( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] )
{
  reductions(V,RTOpPack::ROpNorm1<Scalar>(),norms);
}

template<class Scalar>
inline
void Thyra::norms_2( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] )
{
  reductions(V,RTOpPack::ROpNorm2<Scalar>(),norms);
}

template<class Scalar>
inline
void Thyra::norms_inf( const MultiVectorBase<Scalar>& V, typename Teuchos::ScalarTraits<Scalar>::magnitudeType norms[] )
{
  reductions(V,RTOpPack::ROpNormInf<Scalar>(),norms);
}

#endif // THYRA_MULTI_VECTOR_STD_OPS_DECL_HPP
