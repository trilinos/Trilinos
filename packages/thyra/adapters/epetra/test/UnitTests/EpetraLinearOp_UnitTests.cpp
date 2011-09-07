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


TEUCHOS_UNIT_TEST( EpetraLinearOp, rectangular )
{
  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;

  RCP<Epetra_CrsMatrix> epetraCrsM = getEpetraMatrix(numRows,numCols);

  const RCP<const Thyra::LinearOpBase<double> > epetraOp =
    Thyra::epetraLinearOp(epetraCrsM);

  Thyra::LinearOpTester<double> linearOpTester;

  const bool result = linearOpTester.check(*epetraOp, &out);

  if (!result)
    success = false;
   
}


TEUCHOS_UNIT_TEST( EpetraLinearOp, blocked_op )
{

  using Teuchos::describe;
  using Thyra::block2x2;
  using Thyra::block2x1;
  using Thyra::block1x2;
  using Thyra::LinearOpBase;
  using Thyra::MultiVectorBase;
  using Thyra::createMembers;
  
  // build sub operators
  RCP<const Thyra::LinearOpBase<double> > A00 = Thyra::epetraLinearOp(getEpetraMatrix(4,4,0));
  RCP<const Thyra::LinearOpBase<double> > A01 = Thyra::epetraLinearOp(getEpetraMatrix(4,3,1));
  RCP<const Thyra::LinearOpBase<double> > A02 = Thyra::epetraLinearOp(getEpetraMatrix(4,2,2));
  RCP<const Thyra::LinearOpBase<double> > A10 = Thyra::epetraLinearOp(getEpetraMatrix(3,4,3));
  RCP<const Thyra::LinearOpBase<double> > A11 = Thyra::epetraLinearOp(getEpetraMatrix(3,3,4));
  RCP<const Thyra::LinearOpBase<double> > A12 = Thyra::epetraLinearOp(getEpetraMatrix(3,2,5));
  RCP<const Thyra::LinearOpBase<double> > A20 = Thyra::epetraLinearOp(getEpetraMatrix(2,4,6));
  RCP<const Thyra::LinearOpBase<double> > A21 = Thyra::epetraLinearOp(getEpetraMatrix(2,3,8));
  RCP<const Thyra::LinearOpBase<double> > A22 = Thyra::epetraLinearOp(getEpetraMatrix(2,2,8));
  
  out << "Sub operators built" << std::endl;

  {
     // build composite operator
     RCP<const LinearOpBase<double> > A = block2x2<double>(
       block2x2<double>(A00, A01, A10, A11),   block2x1<double>(A02,A12),
       block1x2<double>(A20, A21),             A22
       );
   
     out << "First composite operator built" << std::endl;
     
     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);
     
     Thyra::randomize(-1.0, 1.0, x.ptr());
   
     out << "A = \n" << describe(*A, Teuchos::VERB_HIGH) << std::endl;
     out << "x = \n" << describe(*x, Teuchos::VERB_HIGH) << std::endl;
     out << "y = \n" << describe(*y, Teuchos::VERB_HIGH) << std::endl;
     
     // perform a matrix vector multiply
     Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr());

     out << "First composite operator completed" << std::endl;
  }

  {
     RCP<const LinearOpBase<double> > A = block2x2<double>(
       A11, block1x2<double>(A10,A12),block2x1<double>(A01,A21),
       block2x2<double>(A00,A02,A20,A22));
     
     out << "Second composite operator built" << std::endl;
     
     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);
     
     Thyra::randomize(-1.0, 1.0, x.ptr());
   
     out << "A = \n" << describe(*A, Teuchos::VERB_HIGH) << std::endl;
     out << "x = \n" << describe(*x, Teuchos::VERB_HIGH) << std::endl;
     out << "y = \n" << describe(*y, Teuchos::VERB_HIGH) << std::endl;
     
     // perform a matrix vector multiply
     Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr());

     out << "Second composite operator completed" << std::endl;
  }

  out << "Test complete" << std::endl;

}


} // namespace
