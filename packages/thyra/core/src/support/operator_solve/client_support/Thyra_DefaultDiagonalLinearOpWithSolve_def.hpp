// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_DefaultDiagonalLinearOpWithSolve_decl.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_TestingTools.hpp" // ToDo: I need to have a better way to get eps()!
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultDiagonalLinearOpWithSolve<Scalar>::DefaultDiagonalLinearOpWithSolve()
{}


template<class Scalar>
DefaultDiagonalLinearOpWithSolve<Scalar>::DefaultDiagonalLinearOpWithSolve(
  const RCP<const VectorBase<Scalar> >   &diag
  )
{
  this->initialize(diag);
}


// protected


// Overridden from LinearOpWithSolveBase


template<class Scalar>
bool
DefaultDiagonalLinearOpWithSolve<Scalar>::solveSupportsImpl(
  EOpTransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return (ST::isComplex ? M_trans==NOTRANS || M_trans==TRANS : true);
}


template<class Scalar>
bool
DefaultDiagonalLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& /* solveMeasureType */) const
{
  return this->solveSupportsImpl(M_trans); // I am a direct solver!
}


template<class Scalar>
SolveStatus<Scalar>
DefaultDiagonalLinearOpWithSolve<Scalar>::solveImpl(
  const EOpTransp transp,
  const MultiVectorBase<Scalar> &B,
  const Ptr<MultiVectorBase<Scalar> > &X,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  ) const
{

#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(this->solveSupportsImpl(transp));
#else
  (void)transp;
#endif

  typedef Teuchos::ScalarTraits<Scalar> ST;

  assign(X, ST::zero());
  
  const Ordinal numCols = B.domain()->dim();
  SolveStatus<Scalar> overallSolveStatus;

  for (Ordinal col_j = 0; col_j < numCols; ++col_j) {

    const RCP<const VectorBase<Scalar> > b = B.col(col_j);
    const RCP<VectorBase<Scalar> > x = X->col(col_j);

    ele_wise_divide( ST::one(), *b, *this->getDiag(), x.ptr() );

  }

  SolveStatus<Scalar> solveStatus;
  solveStatus.solveStatus =
    (nonnull(solveCriteria) && !solveCriteria->solveMeasureType.useDefault()
      ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNKNOWN );
  solveStatus.achievedTol = SolveStatus<Scalar>::unknownTolerance();
  return solveStatus;

}


}	// end namespace Thyra


#endif	// THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_HPP
