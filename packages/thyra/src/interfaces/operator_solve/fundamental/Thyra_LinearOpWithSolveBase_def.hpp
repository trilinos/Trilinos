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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_BASE_DEF_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_BASE_DEF_HPP

#include "Thyra_LinearOpWithSolveBase_decl.hpp"


namespace Thyra {


// Deprecated


template<class Scalar>
bool LinearOpWithSolveBase<Scalar>::solveSupportsConj(EConj conj) const
{
  return solveSupports(applyConjToTrans(conj));
}


template<class Scalar>
bool LinearOpWithSolveBase<Scalar>::solveSupportsSolveMeasureType(
  EConj conj, const SolveMeasureType& solveMeasureType) const
{
  return solveSupportsSolveMeasureType(applyConjToTrans(conj), solveMeasureType);
}


template<class Scalar>
void LinearOpWithSolveBase<Scalar>::solve(
  const EConj conj,
  const MultiVectorBase<Scalar> &B,
  MultiVectorBase<Scalar> *X,
  const int numBlocks,
  const BlockSolveCriteria<Scalar> blockSolveCriteria[],
  SolveStatus<Scalar> blockSolveStatus[]
  ) const
{
  const SolveStatus<Scalar> solveStatus =
    this->solve(applyConjToTrans(conj), B, Teuchos::ptr(X),
      convertBlockSolveCriteriaToSolveCritiera(numBlocks, blockSolveCriteria)
      );
  if (numBlocks) {
    blockSolveStatus[0] = solveStatus;
  }
}


template<class Scalar>
bool LinearOpWithSolveBase<Scalar>::solveTransposeSupportsConj(EConj conj) const
{
  return solveSupports(applyTransposeConjToTrans(conj));
}


template<class Scalar>
bool LinearOpWithSolveBase<Scalar>::solveTransposeSupportsSolveMeasureType(
  EConj conj, const SolveMeasureType& solveMeasureType) const
{
  return solveSupportsSolveMeasureType(applyTransposeConjToTrans(conj),
    solveMeasureType);
}


template<class Scalar>
void LinearOpWithSolveBase<Scalar>::solveTranspose(
  const EConj conj,
  const MultiVectorBase<Scalar> &B,
  MultiVectorBase<Scalar> *X,
  const int numBlocks,
  const BlockSolveCriteria<Scalar> blockSolveCriteria[],
  SolveStatus<Scalar> blockSolveStatus[]
  ) const
{
  const SolveStatus<Scalar> solveStatus =
    this->solve(applyTransposeConjToTrans(conj), B, Teuchos::ptr(X),
      convertBlockSolveCriteriaToSolveCritiera(numBlocks, blockSolveCriteria)
      );
  if (numBlocks) {
    blockSolveStatus[0] = solveStatus;
  }
}


// Protected virtual functions to be overridden by subclasses


template<class Scalar>
bool
LinearOpWithSolveBase<Scalar>::solveSupportsImpl(
  EOpTransp transp) const
{
  return (transp == NOTRANS);
}


template<class Scalar>
bool
LinearOpWithSolveBase<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp transp, const SolveMeasureType& solveMeasureType) const
{
  return (solveSupports(transp) && solveMeasureType.useDefault());
}


// private:


template<class Scalar>
Ptr<const SolveCriteria<Scalar> >
LinearOpWithSolveBase<Scalar>::convertBlockSolveCriteriaToSolveCritiera(
  const int numBlocks,
  const BlockSolveCriteria<Scalar> blockSolveCriteria[]
  )
{
  TEUCHOS_ASSERT(numBlocks == 0 || numBlocks == 1);
  if (numBlocks == 1)
    return Teuchos::ptrFromRef(blockSolveCriteria[0].solveCriteria);
  return Teuchos::null;
}


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_DEF_HPP
