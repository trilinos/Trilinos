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

template <class Scalar>
const typename Teuchos::ScalarTraits<Scalar>::magnitudeType
SolveTolerance<Scalar>::DEFAULT_SOLVE_TOLERANCE = -1;

template <class Scalar>
const typename Teuchos::ScalarTraits<Scalar>::magnitudeType
SolveReturn<Scalar>::UNKNOWN_SOLVE_TOLERANCE = -1;

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveSupports(EConj conj) const
{
  return true;
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveTransposeSupports(EConj conj) const
{
  return true;
}

template <class RangeScalar, class DomainScalar>
SolveReturn<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
LinearOpWithSolveBase<RangeScalar,DomainScalar>::solveTranspose(
  const EConj                           conj
  ,const MultiVectorBase<DomainScalar>  &B
  ,MultiVectorBase<RangeScalar>         *X
  ,const int                            numBlocks
  ,const BlockSolveTolerance<Scalar>    blockSolveTolerances[]
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"LinearOpWithSolveBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::solveTranspose(...): "
    "Error, the concrete subclass described as { " << this->description() << " } "
    " with this->solveTransposeSupports("<<toString(conj)<<")="<<this->solveTransposeSupports(conj)
    << " did not override this function and does not support transposes."
    );
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_HPP
