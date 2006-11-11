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

#ifndef THYRA_INVERSE_LINEAR_OPERATOR_HPP
#define THYRA_INVERSE_LINEAR_OPERATOR_HPP

#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

namespace Thyra {

/** \brief Create an inverse linear operator wrapped as an
 * <tt>DefaultInverseLinearOp</tt> object.
 *
 * \brief  lowsf
 *           [in] The factory used to create the embedded
 *           <tt>LinearOpWithSolveBase</tt> object.
 * \brief  fwdOp
 *           [in] The forward linear operator that is used
 *           to define the <tt>LinearOpWithSolveBase</tt> object.
 *           Warning!  This object will be remembered by the
 *           returned object!  The client can not change
 *           the <tt>*fwdOp.ptr()</tt> object while the returned
 *           inverse operator is still in use.
 * \brief  prevInverseOp
 *           [in] The inverse operator returned from a previous
 *           call to this function.  This object may contain
 *           a great deal of preprocessing data that can be reused
 *           in most cases.
 * 
 * \relates ConstLinearOperator
 */
template<class Scalar>
LinearOperator<Scalar>
inverse(
  LinearOpWithSolveFactoryBase<Scalar>    const&  lowsf
  ,ConstLinearOperator<Scalar>            const&  fwdOp
  ,EThrowOnSolveFailure                           throwOnSolveFailure  = THROW_ON_SOLVE_FAILURE
  ,LinearOperator<Scalar>                 const&  prevInverseOp        = Teuchos::null
  )
{
  return inverse<Scalar>(
    linearOpWithSolve(lowsf,fwdOp.constPtr())
    ,NULL,throwOnSolveFailure
    ,NULL,throwOnSolveFailure
    );
}

} // end namespace Thyra

#endif	// THYRA_INVERSE_LINEAR_OPERATOR_HPP
