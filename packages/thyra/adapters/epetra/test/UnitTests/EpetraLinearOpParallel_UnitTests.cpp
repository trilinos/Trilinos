/*
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
*/

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "EpetraThyraAdaptersTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Unit Tests
//


TEUCHOS_UNIT_TEST( EpetraLinearOpParallel, blocked_op )
{

  using Teuchos::describe;
  using Thyra::block2x2;
  using Thyra::block2x1;
  using Thyra::block1x2;
  using Thyra::LinearOpBase;
  using Thyra::MultiVectorBase;
  using Thyra::createMembers;

  int scale=30;
  
  // build sub operators
  RCP<const Thyra::LinearOpBase<double> > A00 = Thyra::epetraLinearOp(getEpetraMatrix(scale*4,scale*4,0));
  RCP<const Thyra::LinearOpBase<double> > A01 = Thyra::epetraLinearOp(getEpetraMatrix(scale*4,scale*3,1));
  RCP<const Thyra::LinearOpBase<double> > A02 = Thyra::epetraLinearOp(getEpetraMatrix(scale*4,scale*2,2));
  RCP<const Thyra::LinearOpBase<double> > A10 = Thyra::epetraLinearOp(getEpetraMatrix(scale*3,scale*4,3));
  RCP<const Thyra::LinearOpBase<double> > A11 = Thyra::epetraLinearOp(getEpetraMatrix(scale*3,scale*3,4));
  RCP<const Thyra::LinearOpBase<double> > A12 = Thyra::epetraLinearOp(getEpetraMatrix(scale*3,scale*2,5));
  RCP<const Thyra::LinearOpBase<double> > A20 = Thyra::epetraLinearOp(getEpetraMatrix(scale*2,scale*4,6));
  RCP<const Thyra::LinearOpBase<double> > A21 = Thyra::epetraLinearOp(getEpetraMatrix(scale*2,scale*3,8));
  RCP<const Thyra::LinearOpBase<double> > A22 = Thyra::epetraLinearOp(getEpetraMatrix(scale*2,scale*2,8));
  
  out << "Sub operators built" << std::endl;

  typedef double S;
  const RCP<Thyra::PhysicallyBlockedLinearOpBase<S> > M00 = rcp(new Thyra::DefaultBlockedLinearOp<S>());
  M00->beginBlockFill(2,2);
  M00->setBlock(0,0,A00);
  M00->setBlock(0,1,A02);
  M00->setBlock(1,0,A20);
  M00->setBlock(1,1,A22);
  M00->endBlockFill();

  out << "Built M00" << std::endl;

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<S> > M10 = rcp(new Thyra::DefaultBlockedLinearOp<S>());
  M10->beginBlockFill(2,1);
  M10->setBlock(0,0,A01);
  M10->setBlock(1,0,A21);
  M10->endBlockFill();

  out << "Built M10" << std::endl;

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<S> > M01 = rcp(new Thyra::DefaultBlockedLinearOp<S>());
  M01->beginBlockFill(1,2);
  M01->setBlock(0,0,A10);
  M01->setBlock(0,1,A12);
  M01->endBlockFill();

  out << "Built M01" << std::endl;

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<S> > M = rcp(new Thyra::DefaultBlockedLinearOp<S>());
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


} // namespace
