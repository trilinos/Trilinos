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

#ifndef THYRA_EUCLIDEAN_LINEAR_OP_DECL_HPP
#define THYRA_EUCLIDEAN_LINEAR_OP_DECL_HPP

#include "Thyra_OperatorVectorAdapterSupportTypes.hpp"
#include "Thyra_LinearOpDefaultBaseDecl.hpp"

namespace Thyra {

/** \brief Base interface for Euclidean linear operators.
 *
 * Most generic subclass implementations should derive from this base
 * interface since it allows an application-specific definition of the
 * scalar product.  Note that almost every concrete implementation of
 * <tt>MultiVectorBase</tt> should derive from this interface.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class RangeScalar, class DomainScalar>
class EuclideanLinearOpBase : virtual public LinearOpDefaultBase<RangeScalar,DomainScalar> {
public:

  /** @name Pure virtual functions to override in subclasses */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<RangeScalar> > rangeScalarProdVecSpc() const = 0;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<DomainScalar> > domainScalarProdVecSpc() const = 0;

  /** \brief Apply the linear operator to a multi-vector with respect
   * to a Euclidean vector space where the scalar product is the dot
   * product.
   *
   * Preconditions:<ul>
   * <li><tt>this->applySupports(conj)==true</tt>
   * </ul>
   */
  virtual void euclideanApply(
    const EConj                            conj
    ,const MultiVectorBase<DomainScalar>   &X
    ,MultiVectorBase<RangeScalar>          *Y
    ,const RangeScalar                     alpha
    ,const RangeScalar                     beta
    ) const = 0;

  //@}

  /** @name Virtual functions with default implementations. */
  //@{

  /** \brief Apply the linear operator to a multi-vector with respect
   * to a Euclidean vector space where the scalar product is the dot
   * product.
   *
   * Preconditions:<ul>
   * <li><tt>this->applyTransposeSupports(conj)==true</tt>
   * </ul>
   *
   * The default implementation throws an exception with a very good error
   * message.
   */
  virtual void euclideanApplyTranspose(
    const EConj                            conj
    ,const MultiVectorBase<RangeScalar>    &X
    ,MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                    alpha
    ,const DomainScalar                    beta
    ) const;

  //@}

  /** @name Overridden functions from OpBase */
  //@{
  /// Returns <tt>this->rangeScalarProdVecSpc()</tt>
  Teuchos::RefCountPtr<const VectorSpaceBase<RangeScalar> > range() const;
  /// Returns <tt>this->domainScalarProdVecSpc()</tt>
  Teuchos::RefCountPtr<const VectorSpaceBase<DomainScalar> > domain() const;
  //@}

  /** @name Overridden functions from LinearOpBase */
  //@{

  /** \brief Apply the non-transposed linear operator to a multi-vector using
   * an application-specific definition of the scalar product.
   *
   * ToDo: Finish Documentation!
   */
  void apply(
    const EConj                            conj
    ,const MultiVectorBase<DomainScalar>   &X
    ,MultiVectorBase<RangeScalar>          *Y
    ,const RangeScalar                     alpha
    ,const RangeScalar                     beta
    ) const;

  /** \brief Apply the transposed linear operator to a multi-vector using
   * an application-specific definition of the scalar product.
   *
   * ToDo: Finish Documentation!
   */
  void applyTranspose(
    const EConj                            conj
    ,const MultiVectorBase<RangeScalar>    &X
    ,MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                    alpha
    ,const DomainScalar                    beta
    ) const;

  //@}

protected:
  
  void euclidean_apply_impl(
    const EConj                            conj
    ,const MultiVectorBase<DomainScalar>   &X
    ,MultiVectorBase<RangeScalar>          *Y
    ,const RangeScalar                     alpha
    ,const RangeScalar                     beta
    ) const;
  
  void euclidean_applyTranspose_impl(
    const EConj                            conj
    ,const MultiVectorBase<RangeScalar>    &X
    ,MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                    alpha
    ,const DomainScalar                    beta
    ) const;

}; // end class EuclideanLinearOpBase<Scalar>

/** \brief Call <tt>EuclideanLinearOpBase<Scalar>::euclideanApply()</tt> as a
 *    global function call (for a single scalar type).
 *
 * Calls <tt>M.euclideanApply(...,X,Y,alpha,beta)</tt> or
 * <tt>M.euclideanApplyTranspose(...,X,Y,alpha,beta)</tt>.
 *
 * \relates EuclideanLinearOpBase
 */
template<class Scalar>
inline void euclideanApply(
  const EuclideanLinearOpBase<Scalar>        &M
  ,const ETransp                             M_trans
  ,const MultiVectorBase<Scalar>             &X
  ,MultiVectorBase<Scalar>                   *Y
  ,const Scalar                              alpha
  ,const Scalar                              beta
  )
{
  if(real_trans(M_trans)==NOTRANS) {
    M.euclideanApply(transToConj(M_trans),X,Y,alpha,beta);
  }
  else {
    M.euclideanApplyTranspose(transToConj(M_trans),X,Y,alpha,beta);
  }
}

} // namespace Thyra

#endif // THYRA_EUCLIDEAN_LINEAR_OP_DECL_HPP
