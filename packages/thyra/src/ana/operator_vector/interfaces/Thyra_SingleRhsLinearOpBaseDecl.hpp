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

#ifndef THYRA_SINGLE_RHS_LINEAR_OP_BASE_DECL_HPP
#define THYRA_SINGLE_RHS_LINEAR_OP_BASE_DECL_HPP

#include "Thyra_OpBaseDecl.hpp"

namespace Thyra {

/** \brief Base class for linear operators that can only implement a single
 * RHS vector apply operation.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
class SingleRhsLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

protected:

  /** @name Pure virtual functions (must be overridden by subclass) */
  //@{

  /** \brief Apply the linear operator (or its transpose) to a vector:
   * <tt>y = alpha*op(M)*x + beta*y</tt>.
   *
   * @param  M_trans
   *                [in] Determines whether the transposed or non-transposed
   *                operator is applied as:
   *                <ul>
   *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
   *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
   *                </ul>
   *                where <tt>M == *this</tt>
   * @param  x      [in] The right hand side vector 
   * @param  y      [in/out] The target vector being transformed
   * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
   * @param  beta   [in] The multiplier for the target vector <tt>y</tt>.
   *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
   * 
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>this->opSupported(M_trans)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
   * <li> <tt>y->space()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>x.space()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>y</tt> can not alias <tt>x</tt>.  It is up to the client to ensure that <tt>y</tt>
   *      and <tt>x</tt> are distinct since in general this can not be verified by the implementation until,
   *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
   * </ul>
   *
   * Postconditions:<ul>
   * <li> Is it not obvious?  After the function returns the vector <tt>y</tt>
   *      is transformed as indicated above.
   * </ul>
   */
  virtual void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const = 0;

  //@}

};	// end class LinearOpBase

}	// end namespace Thyra

#endif	// THYRA_SINGLE_RHS_LINEAR_OP_BASE_DECL_HPP
