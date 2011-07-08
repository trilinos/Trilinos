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

#ifndef THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"


namespace Thyra {


/** \brief Base interface for linear operators with a solve that can be
 * accessed as sub-blocks.
 *
 * \ingroup Thyra_Op_Solve_extended_interfaces_code_grp
 *
 * ToDo: Finish Documentation.
 */
template<class Scalar>
class BlockedLinearOpWithSolveBase
  : virtual public LinearOpWithSolveBase<Scalar>,
    virtual public BlockedLinearOpBase<Scalar>
  
{
public:

  /** \brief . */
  virtual Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLOWSBlock(const int i, const int j) = 0; 

  /** \brief . */
  virtual Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
  getLOWSBlock(const int i, const int j) const = 0; 

private:
  
  // Not defined and not to be called
  BlockedLinearOpWithSolveBase<Scalar>&
  operator=(const BlockedLinearOpWithSolveBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
