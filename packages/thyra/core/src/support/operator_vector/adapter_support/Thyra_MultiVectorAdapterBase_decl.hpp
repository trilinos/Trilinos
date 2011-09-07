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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_ADAPTER_BASE_DECL_HPP
#define THYRA_MULTI_VECTOR_ADAPTER_BASE_DECL_HPP

#include "Thyra_MultiVectorDefaultBase.hpp"


namespace Thyra {


/** \brief Forward decl. */
template<class Scalar> class ScalarProdVectorSpaceBase;


/** \brief Node subclass for MultiVectorBase subclasses that allows the
 * insertion of an application defined scalar product.
 *
 * Most concrete MultiVector adapter subclasses should derive from this base
 * subclass in order to allow for the incorporate of application-defined
 * scalar products.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class MultiVectorAdapterBase : virtual public MultiVectorDefaultBase<Scalar>
{
public:

  /** @name Pure virtual functions to override in subclasses */
  //@{

  /** \brief . */
  virtual RCP<const ScalarProdVectorSpaceBase<Scalar> >
  rangeScalarProdVecSpc() const = 0;

  /** \brief . */
  virtual RCP<const ScalarProdVectorSpaceBase<Scalar> >
  domainScalarProdVecSpc() const = 0;

  /** \brief Apply the linear operator to a multi-vector with respect
   * to a Euclidean vector space where the scalar product is the dot
   * product.
   *
   * Preconditions:<ul>
   * <li><tt>this->applySupports(conj)==true</tt>
   * </ul>
   */
  virtual void euclideanApply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar>   &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const = 0;

  //@}

  /** @name Overridden functions from LinearOp */
  //@{
  /// Returns <tt>this->rangeScalarProdVecSpc()</tt>
  RCP<const VectorSpaceBase<Scalar> > range() const;
  /// Returns <tt>this->domainScalarProdVecSpc()</tt>
  RCP<const VectorSpaceBase<Scalar> > domain() const;
  //@}

protected:

  /** @name Overridden protected functions from LinearOpBase */
  //@{
  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief .  */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;
  //@}

};


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_ADAPTER_BASE_DECL_HPP
