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

#include "tPCDStrategy_tpetra.hpp"

#include <string>

// Teko-Package includes
#include "Teko_Utilities.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

// Teko includes
#include "Teko_PCDStrategy.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StaticRequestCallback.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::tpetraLinearOp;

RCP<Thyra::LinearOpBase<ST> > buildExampleOp(int type,RCP<const Teuchos::Comm<int> > comm,ST scale)
{
   const RCP<const Tpetra::Map<LO,GO,NT> > map = rcp(new const Tpetra::Map<LO,GO,NT>(3,0,comm));
   const RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > mat = Tpetra::createCrsMatrix<ST,LO,GO,NT> (map,3);

   ST values[3];
   GO indices[3] = { 0, 1, 2 };

   switch(type) {
   case 3:
      values[0] = 9.0; values[1] = 8.0; values[2] = 7.0;
      mat->insertGlobalValues(0,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));

      values[0] = 1.0; values[1] = 2.0; values[2] = 3.0;
      mat->insertGlobalValues(1,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));

      values[0] = 4.0; values[1] = 5.0; values[2] = 6.0;
      mat->insertGlobalValues(2,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
      break;
   case 2:
      values[0] = 5.0; values[1] = 0.0; values[2] = 6.0;
      mat->insertGlobalValues(0,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));

      values[0] = 2.0; values[1] = 6.0; values[2] = 0.0;
      mat->insertGlobalValues(1,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));

      values[0] = 0.0; values[1] = 3.0; values[2] = 7.0;
      mat->insertGlobalValues(2,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
      break;
   case 1:
      values[0] = 1.0; values[1] = 0.0; values[2] = 2.0;
      mat->insertGlobalValues(0,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));

      values[0] = 2.0; values[1] = 5.0; values[2] = 0.0;
      mat->insertGlobalValues(1,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));

      values[0] = 0.0; values[1] = 3.0; values[2] = 1.0;
      mat->insertGlobalValues(2,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
   default:
      break;
   }

   mat->scale(scale);
   mat->fillComplete();

   return Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(mat->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(mat->getDomainMap()),mat);
}

void tPCDStrategy_tpetra::initializeTest()
{
   tolerance_ = 1.0e-15;
}

int tPCDStrategy_tpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tPCDStrategy_tpetra";

   status = test_PCDStrategy(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"PCDStrategy\" ... PASSED","   \"PCDStrategy\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tPCDStrategy_tpetra...PASSED","tPCDStrategy_tpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tPCDStrategy_tpetra...FAILED");
   }

   return failcount;
}

bool tPCDStrategy_tpetra::test_PCDStrategy(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(tolerance_);

   Teko::LinearOp M1 = buildExampleOp(1,GetComm_tpetra(),1.0);
   Teko::LinearOp M2 = buildExampleOp(2,GetComm_tpetra(),1.0);
   Teko::LinearOp A = Thyra::nonconstBlock2x2<ST>(buildExampleOp(2,GetComm_tpetra(),1.0),
                                      buildExampleOp(1,GetComm_tpetra(),2.0),
                                      buildExampleOp(1,GetComm_tpetra(),3.0),
                                      buildExampleOp(2,GetComm_tpetra(),4.0));
   Teko::BlockedLinearOp blkA = Teko::toBlockedLinearOp(A);

   Teko::LinearOp laplace = Teko::explicitMultiply(M2,M1);
   Teko::LinearOp Fp = M1;
   Teko::LinearOp Qp = buildExampleOp(3,GetComm_tpetra(),1.0);

   std::string pcdStr      = Teko::NS::PCDStrategy::getPCDString();
   std::string presLapStr  = Teko::NS::PCDStrategy::getPressureLaplaceString();
   std::string presMassStr = Teko::NS::PCDStrategy::getPressureMassString();

   RCP<Teko::StaticRequestCallback<Teko::LinearOp> > pcdCb
         = rcp(new Teko::StaticRequestCallback<Teko::LinearOp>(pcdStr,Fp));
   RCP<Teko::StaticRequestCallback<Teko::LinearOp> > lapCb
         = rcp(new Teko::StaticRequestCallback<Teko::LinearOp>(presLapStr,laplace));
   RCP<Teko::StaticRequestCallback<Teko::LinearOp> > massCb
         = rcp(new Teko::StaticRequestCallback<Teko::LinearOp>(presMassStr,Qp));

   RCP<Teko::RequestHandler> rh = rcp(new Teko::RequestHandler);
   rh->addRequestCallback(pcdCb);
   rh->addRequestCallback(lapCb);
   rh->addRequestCallback(massCb);

   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   RCP<Teko::InverseFactory> inverseFact = invLib->getInverseFactory("Ifpack2");

   // test the inverse of the 00 block
   {
      bool result;
      std::stringstream ss;
      Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss),"      |||");

      Teko::LinearOp invA00 = inverseFact->buildInverse(M2);
      Teko::BlockPreconditionerState state;
      Teko::NS::PCDStrategy strategy(inverseFact,inverseFact);
      strategy.setRequestHandler(rh);

      // build the state object
      // state.addLinearOp(presLapStr,  laplace);
      // state.addLinearOp(pcdStr,      Fp);
      // state.addLinearOp(presMassStr, Qp);

      // test the 0,0 block inverses: hat
      Teko::LinearOp hatInvA00 = strategy.getHatInvA00(blkA,state);
      ss.str("");
      result = tester.compare( *hatInvA00, *invA00, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tPCDStrategy_tpetra::test_PCDStrategy " << toString(status)
                      << " : Comparing invHatA00 operators");
      if(not result || verbosity>=10)
         os << ss.str();

      // test the 0,0 block inverses: tilde
      Teko::LinearOp tildeInvA00 = strategy.getTildeInvA00(blkA,state);
      ss.str("");
      result = tester.compare( *tildeInvA00, *invA00, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tPCDStrategy_tpetra::test_PCDStrategy " << toString(status)
                      << " : Comparing invTildeA00 operators");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   // test the inverse of the Schur complement
   {
      bool result;
      std::stringstream ss;
      Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss),"      |||");

      Teko::LinearOp invLaplace = inverseFact->buildInverse(laplace);
      Teko::LinearOp iQp = Teko::getInvDiagonalOp(Qp);
      Teko::LinearOp invSchur = multiply(iQp,Fp,invLaplace);

      Teko::BlockPreconditionerState state;
      Teko::NS::PCDStrategy strategy(inverseFact,inverseFact);
      strategy.setRequestHandler(rh);

      // build the state object
      // state.addLinearOp(presLapStr,  laplace);
      // state.addLinearOp(pcdStr,      Fp);
      // state.addLinearOp(presMassStr, Qp);

      // test the 0,0 block inverses: hat
      Teko::LinearOp invS_strategy = strategy.getInvS(blkA,state);
      ss.str("");
      result = tester.compare( *invS_strategy, *invSchur, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tPCDStrategy_tpetra::test_PCDStrategy " << toString(status)
                      << " : Comparing invS operators");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
