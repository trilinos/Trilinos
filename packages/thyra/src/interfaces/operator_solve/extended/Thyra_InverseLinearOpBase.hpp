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

#ifndef THYRA_INVERSE_LINEAR_OP_BASE_HPP
#define THYRA_INVERSE_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"


namespace Thyra {


/** \brief Base interface for <ttLinearOpBase</tt> objects that are implemented
 * in terms of the solve function on a <tt>LinearOpWithSolveBase</tt> object.
 *
 * \ingroup Thyra_Op_Solve_extended_interfaces_code_grp
 */
template<class Scalar>
class InverseLinearOpBase : virtual public LinearOpBase<Scalar>
{
public:

  /** \brief Determine if the underlying <tt>LinearOpWithSolveBase</tt> is
   * const-only or not.
   */
  virtual bool isLowsConst() const = 0;

  /** \brief Extra a non-const view of the underlying
   * <tt>LinearOpWithSolveBase</tt> object.
   */
  virtual Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLows() = 0; 

  /** \brief Extra a const view of the underlying
   * <tt>LinearOpWithSolveBase</tt> object.
   */
  virtual Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
  getLows() const = 0; 

};


} // namespace Thyra


#endif // THYRA_INVERSE_LINEAR_OP_BASE_HPP
