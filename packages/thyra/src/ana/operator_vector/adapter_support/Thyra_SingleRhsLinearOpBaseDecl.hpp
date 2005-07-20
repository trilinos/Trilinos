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

#include "Thyra_SingleScalarLinearOpBaseDecl.hpp"

namespace Thyra {

/** \brief Base class for linear operators that can only implement a single
 * RHS vector apply operation and only support one scalar type.
 *
 * This class is meant to provide an easier way for subclasses to provide
 * implementations for the multi-vector version of
 * <tt>SingleScalarLinearOpBase::apply()</tt> and is not meant to be used as
 * an client interface.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
class SingleRhsLinearOpBase : virtual public SingleScalarLinearOpBase<Scalar> {
public:

  /** \brief . */
  using SingleScalarLinearOpBase<Scalar>::apply;

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

  /** \brief Apply the linear operator (or its transpose) to single vector
   * arguments.
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
