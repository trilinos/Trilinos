// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_BASE_DEF_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_BASE_DEF_HPP

#include "Thyra_LinearOpWithSolveBase_decl.hpp"


namespace Thyra {


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


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_DEF_HPP
