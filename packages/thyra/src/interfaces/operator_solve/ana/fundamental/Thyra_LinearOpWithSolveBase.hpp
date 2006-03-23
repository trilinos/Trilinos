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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace Thyra {

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveSupportsConj(EConj conj) const
{
  return true;
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveTransposeSupportsConj(EConj conj) const
{
  return false;
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveSupportsSolveMeasureType(EConj conj, const SolveMeasureType& solveMeasureType) const
{
  return solveMeasureType.useDefault();
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveTransposeSupportsSolveMeasureType(EConj conj, const SolveMeasureType& solveMeasureType) const
{
  return solveMeasureType.useDefault();
}

template <class RangeScalar, class DomainScalar>
void LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveTranspose(
  const EConj                                   conj
  ,const MultiVectorBase<DomainScalar>          &B
  ,MultiVectorBase<RangeScalar>                 *X
  ,const int                                    numBlocks
  ,const BlockSolveCriteria<PromotedScalar>     blockSolveCriteria[]
  ,SolveStatus<PromotedScalar>                  blockSolveStatus[]
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"LinearOpWithSolveBase<"<<Teuchos::ScalarTraits<RangeScalar>::name()<<","<<Teuchos::ScalarTraits<DomainScalar>::name()<<">::solveTranspose(...): "
    "Error, the concrete subclass described as { " << this->description() << " } "
    " with this->solveTransposeSupportsConj("<<toString(conj)<<")="<<this->solveTransposeSupportsConj(conj)
    << " did not override this function and does not support transposes."
    );
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_HPP
