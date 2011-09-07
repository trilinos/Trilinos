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
