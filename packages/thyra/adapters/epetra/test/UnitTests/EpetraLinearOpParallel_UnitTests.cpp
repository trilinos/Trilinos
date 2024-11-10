// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "EpetraLinearOpTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Thyra {


//
// Unit Tests
//


TEUCHOS_UNIT_TEST( EpetraLinearOpParallel, blocked_op )
{

  using Teuchos::describe;

  const int scale=30;
  
  // build sub operators
  RCP<const LinearOpBase<double> > A00 = epetraLinearOp(getEpetraMatrix(scale*4,scale*4,0));
  RCP<const LinearOpBase<double> > A01 = epetraLinearOp(getEpetraMatrix(scale*4,scale*3,1));
  RCP<const LinearOpBase<double> > A02 = epetraLinearOp(getEpetraMatrix(scale*4,scale*2,2));
  RCP<const LinearOpBase<double> > A10 = epetraLinearOp(getEpetraMatrix(scale*3,scale*4,3));
  RCP<const LinearOpBase<double> > A11 = epetraLinearOp(getEpetraMatrix(scale*3,scale*3,4));
  RCP<const LinearOpBase<double> > A12 = epetraLinearOp(getEpetraMatrix(scale*3,scale*2,5));
  RCP<const LinearOpBase<double> > A20 = epetraLinearOp(getEpetraMatrix(scale*2,scale*4,6));
  RCP<const LinearOpBase<double> > A21 = epetraLinearOp(getEpetraMatrix(scale*2,scale*3,8));
  RCP<const LinearOpBase<double> > A22 = epetraLinearOp(getEpetraMatrix(scale*2,scale*2,8));
  
  out << "Sub operators built" << std::endl;

  typedef double S;
  const RCP<PhysicallyBlockedLinearOpBase<S> > M00 = defaultBlockedLinearOp<S>();
  M00->beginBlockFill(2,2);
  M00->setBlock(0,0,A00);
  M00->setBlock(0,1,A02);
  M00->setBlock(1,0,A20);
  M00->setBlock(1,1,A22);
  M00->endBlockFill();

  out << "Built M00" << std::endl;

  const RCP<PhysicallyBlockedLinearOpBase<S> > M10 = defaultBlockedLinearOp<S>();
  M10->beginBlockFill(2,1);
  M10->setBlock(0,0,A01);
  M10->setBlock(1,0,A21);
  M10->endBlockFill();

  out << "Built M10" << std::endl;

  const RCP<PhysicallyBlockedLinearOpBase<S> > M01 = defaultBlockedLinearOp<S>();
  M01->beginBlockFill(1,2);
  M01->setBlock(0,0,A10);
  M01->setBlock(0,1,A12);
  M01->endBlockFill();

  out << "Built M01" << std::endl;

  const RCP<PhysicallyBlockedLinearOpBase<S> > M = defaultBlockedLinearOp<S>();
  M->beginBlockFill(2,2);
  M->setBlock(0,0,A11);
  M->setBlock(0,1,M01);
  M->setBlock(1,0,M10);
  M->setBlock(1,1,M00);
  M->endBlockFill();

  out << "Built M" << std::endl;

  out << "M = " << describe(*M, Teuchos::VERB_MEDIUM);
  out << "M->range() = " << describe(*M->range(), Teuchos::VERB_MEDIUM);
  out << "M->domain() = " << describe(*M->range(), Teuchos::VERB_MEDIUM);

  out << "Test complete" << std::endl;

}


} // namespace Thyra
