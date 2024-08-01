// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_MultiContinuation_MultiVecConstraint.H"

void
LOCA::BorderedSolver::AbstractStrategy::setMatrixBlocksMultiVecConstraint(
          const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockB,
      const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  Teuchos::RCP<LOCA::MultiContinuation::MultiVecConstraint> con =
    Teuchos::rcp(new LOCA::MultiContinuation::MultiVecConstraint(blockB));

  setMatrixBlocks(op, blockA, con, blockC);
}
