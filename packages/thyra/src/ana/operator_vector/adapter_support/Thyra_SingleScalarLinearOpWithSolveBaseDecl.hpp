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

#ifndef THYRA_SINGLE_SCALAR_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
#define THYRA_SINGLE_SCALAR_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"
#include "Thyra_SingleScalarLinearOpBaseDecl.hpp"

namespace Thyra {

/** \brief . */
template <class Scalar>
class SingleScalarLinearOpWithSolveBase
  : virtual public LinearOpWithSolveBase<Scalar>
  , virtual protected SingleScalarLinearOpBase<Scalar>
{
public:
  
  /** @name Overridden from LinearOpWithSolveBase */
  //@{

  /** \brief . */
  bool solveSupports(EConj conj) const;

  /** \brief . */
  bool solveTransposeSupports(EConj conj) const;

  /** \brief . */
  SolveReturn<Scalar> solve(
    const EConj                           conj
    ,const MultiVectorBase<Scalar>        &B
    ,MultiVectorBase<Scalar>              *X
    ,const int                            numBlocks
    ,const BlockSolveTolerance<Scalar>    blockSolveTolerances[]
    ) const;

  /** \brief . */
  SolveReturn<Scalar> solveTranspose(
    const EConj                           conj
    ,const MultiVectorBase<Scalar>        &B
    ,MultiVectorBase<Scalar>              *X
    ,const int                            numBlocks
    ,const BlockSolveTolerance<Scalar>    blockSolveTolerances[]
    ) const;

  //@}

protected:

  /** @name Protected pure virtual functions to be overridden by subclasses. */
  //@{

  /** \brief . */
  virtual bool solveSupported(ETransp M_trans) const = 0;

  /** \brief . */
  virtual SolveReturn<Scalar> solve(
    const ETransp                         M_trans
    ,const MultiVectorBase<Scalar>        &B
    ,MultiVectorBase<Scalar>              *X
    ,const int                            numBlocks
    ,const BlockSolveTolerance<Scalar>    blockSolveTolerances[]
    ) const = 0;

  //@}

};

} // namespace Thyra

#endif // THYRA_SINGLE_SCALAR_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
