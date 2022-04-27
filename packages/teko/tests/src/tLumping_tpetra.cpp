/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include "tLumping_tpetra.hpp"

#include <string>

// Tpetra includes

// Teko-Package includes
#include "Teko_Utilities.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::tpetraLinearOp;

void tLumping_tpetra::initializeTest()
{
   tolerance_ = 1.0e-15;
}

int tLumping_tpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLumping_tpetra";

   status = test_lumping(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"lumping\" ... PASSED","   \"lumping\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_invLumping(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"invLumping\" ... PASSED","   \"invLumping\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tLumping_tpetra...PASSED","tLumping_tpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tLumping_tpetra...FAILED");
   }

   return failcount;
}

bool tLumping_tpetra::test_lumping(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(tolerance_);

   Tpetra::Map<LO,GO,NT> map(100,0,GetComm_tpetra());

   // A matrix...to be lumped
   Tpetra::CrsMatrix<ST,LO,GO,NT> A(rcpFromRef(map),5);
   GO indices[5];
   ST values[5] = {1,2,3,4,5};
   for(size_t i=0;i<A.getLocalNumRows()-5;i++) {
      GO index = A.getRowMap()->getGlobalElement(i);
      for(GO j=0;j<5;j++)
         indices[j] = A.getRowMap()->getGlobalElement(i+j);
      A.insertGlobalValues(index,Teuchos::ArrayView<GO>(indices,5),Teuchos::ArrayView<ST>(values,5));
   }
   for(size_t i=A.getLocalNumRows()-5;i<A.getLocalNumRows();i++) {
      GO index = A.getRowMap()->getGlobalElement(i);
      for(LO j=0;j<5;j++)
         indices[j] = A.getRowMap()->getGlobalElement(j);
      A.insertGlobalValues(index,Teuchos::ArrayView<GO>(indices,5),Teuchos::ArrayView<ST>(values,5));
   }
   A.fillComplete();

   // B matrix...already lumped
   Tpetra::CrsMatrix<ST,LO,GO,NT> B(rcpFromRef(map),1);
   ST number[1] = {15.0};
   for(size_t i=0;i<B.getLocalNumRows();i++) {
      GO index[1] = {B.getRowMap()->getGlobalElement(i)};
      B.insertGlobalValues(index[0],Teuchos::ArrayView<GO>(index,1),Teuchos::ArrayView<ST>(number,1));
   }
   B.fillComplete();

   Teko::LinearOp pA = Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(A.getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(A.getDomainMap()),rcpFromRef(A));
   Teko::LinearOp pB = Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(B.getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(B.getDomainMap()),rcpFromRef(B));
   Teko::LinearOp lumpedA = getLumpedMatrix(pA);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *pB, *lumpedA, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tLummping::test_lumping "
             << ": Testing basic lumping functionality");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   return allPassed;
}

bool tLumping_tpetra::test_invLumping(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(tolerance_);

   Tpetra::Map<LO,GO,NT> map(100,0,GetComm_tpetra());

   // A matrix...to be lumped
   Tpetra::CrsMatrix<ST,LO,GO,NT> A(rcpFromRef(map),5);
   GO indices[5];
   ST values[5] = {1,2,3,4,5};
   for(size_t i=0;i<A.getLocalNumRows()-5;i++) {
      GO index = A.getRowMap()->getGlobalElement(i);
      for(GO j=0;j<5;j++)
         indices[j] = A.getRowMap()->getGlobalElement(i+j);
      A.insertGlobalValues(index,Teuchos::ArrayView<GO>(indices,5),Teuchos::ArrayView<ST>(values,5));
   }
   for(size_t i=A.getLocalNumRows()-5;i<A.getLocalNumRows();i++) {
      GO index = A.getRowMap()->getGlobalElement(i);
      for(LO j=0;j<5;j++)
         indices[j] = A.getRowMap()->getGlobalElement(j);
      A.insertGlobalValues(index,Teuchos::ArrayView<GO>(indices,5),Teuchos::ArrayView<ST>(values,5));
   }
   A.fillComplete();

   // B matrix...already lumped
   Tpetra::CrsMatrix<ST,LO,GO,NT> B(rcpFromRef(map),1);
   ST number[1] = {1.0/15.0};
   for(size_t i=0;i<B.getLocalNumRows();i++) {
      GO index[1] = {B.getRowMap()->getGlobalElement(i)};
      B.insertGlobalValues(index[0],Teuchos::ArrayView<GO>(index,1),Teuchos::ArrayView<ST>(number,1));
   }
   B.fillComplete();

   Teko::LinearOp pA = Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(A.getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(A.getDomainMap()),rcpFromRef(A));
   Teko::LinearOp pB = Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(B.getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(B.getDomainMap()),rcpFromRef(B));
   Teko::LinearOp lumpedA = getInvLumpedMatrix(pA);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *pB, *lumpedA, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tLummping::test_invLumping "
             << ": Testing basic inv lumping functionality");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
