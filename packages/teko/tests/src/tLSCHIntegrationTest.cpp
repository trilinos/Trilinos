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

#include "tLSCHIntegrationTest.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// ML includes
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

// Teko-Package includes
#include "Teko_EpetraHelpers.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_Utilities.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Test-rig
#include "Test_Utils.hpp"

#include "Teko_InverseFactory.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::epetraLinearOp;

void tLSCHIntegrationTest::initializeTest()
{
   tolerance_ = 1.0e-7;
}

int tLSCHIntegrationTest::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLSCHIntegrationTest";

   status = test_hScaling(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"hScaling\" ... PASSED","   \"hScaling\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tLSCHIntegrationTest...PASSED","tLSCHIntegrationTest...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tLSCHIntegrationTest...FAILED");
   }

   return failcount;
}

bool tLSCHIntegrationTest::test_hScaling(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;
   
   RCP<const Epetra_Comm> comm = GetComm();

   // build some operators
   Teko::LinearOp F = Teko::Test::build2x2(*comm,1,2,2,1);
   Teko::LinearOp G = Teko::Test::build2x2(*comm,1,-1,-3,1);
   Teko::LinearOp D = Teko::Test::build2x2(*comm,1,-3,-1,1);

   double diag[2];
 
   diag[0] = 1.0/3.0; diag[1] = 1.0/2.0;
   Teko::LinearOp M = Teko::Test::DiagMatrix(2,diag,"M");

   diag[0] = 5.0; diag[1] = 9.0;
   Teko::LinearOp H = Teko::Test::DiagMatrix(2,diag,"H");

   Teko::LinearOp A = Thyra::block2x2<double>(F,G,D,Teuchos::null);

   Teko::LinearOp exact;
   {
      // build some operators
      Teko::LinearOp D0 = Teko::Test::build2x2(*comm,-1.0/3.0,2.0/3.0,2.0/3.0,-1.0/3.0);
      Teko::LinearOp D1 = Teko::Test::build2x2(*comm,-1.5,-3.0,-3.0,-5.5);
      Teko::LinearOp U  = Teko::Test::build2x2(*comm,-0.5,-1.5,-0.5,-0.5);
      
      exact = Thyra::block2x2<double>(D0,U,Teuchos::null,D1);
   }

   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");
   RCP<Teko::NS::InvLSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(invFact,M));
   strategy->setHScaling(Teko::getDiagonal(H));
   strategy->setUseFullLDU(false);

   RCP<Teko::BlockPreconditionerFactory> precFact = rcp(new Teko::NS::LSCPreconditionerFactory(strategy));
   RCP<Teko::BlockPreconditionerState> bps = Teuchos::rcp_dynamic_cast<Teko::BlockPreconditionerState>(precFact->buildPreconditionerState());
   Teko::LinearOp prec = precFact->buildPreconditionerOperator(A,*bps);

   Teko::BlockedLinearOp bA = Teko::toBlockedLinearOp(A);
   std::stringstream ss;
   ss << "invF = " << Teuchos::describe(*strategy->getInvF(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "invBQBt = " << Teuchos::describe(*strategy->getInvBQBt(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "invF = " << Teuchos::describe(*strategy->getInvBHBt(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "invMass = " << Teuchos::describe(*strategy->getInvMass(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "HScaling = " << Teuchos::describe(*strategy->getHScaling(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "prec = " << Teuchos::describe(*prec,Teuchos::VERB_EXTREME) << std::endl;

   // construct a couple of vectors
   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*GetComm()));
   Epetra_Vector ea(*map),eb(*map);
   const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea,eb,prec->domain());
   const RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers(prec->range(),1); 
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] = 0.0;
   Thyra::apply(*prec,Thyra::NOTRANS,*x,y.ptr());
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 0.0; eb[1] = 0.0;
   Thyra::apply(*prec,Thyra::NOTRANS,*x,y.ptr());
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
   ea[0] = 0.0; ea[1] = 0.0; eb[0] = 1.0; eb[1] = 0.0;
   Thyra::apply(*prec,Thyra::NOTRANS,*x,y.ptr());
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
   ea[0] = 0.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] = 1.0;
   Thyra::apply(*prec,Thyra::NOTRANS,*x,y.ptr());
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
  
   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss),"|||");
   const bool result = tester.compare( *prec, *exact, Teuchos::ptrFromRef(fos) );
   TEST_ASSERT(result,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing preconditioner to exactly computed version");

   const bool result2 = tester.compare( *H, *strategy->getHScaling(bA,*bps), Teuchos::ptrFromRef(fos) );
   TEST_ASSERT(result2,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing scaling of H operator");

   const bool result3 = tester.compare( *build2x2(*comm,0.208333333333333,0.375, 0.375, 0.875),
     *strategy->getInvBQBt(bA,*bps), Teuchos::ptrFromRef(fos) );
   TEST_ASSERT(result3,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing inv(BQBt) operator");

   const bool result4 = tester.compare( *build2x2(*comm, 0.077777777777778, 0.177777777777778, 0.177777777777778, 0.477777777777778),
     *strategy->getInvBHBt(bA,*bps), Teuchos::ptrFromRef(fos) );
   TEST_ASSERT(result4,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing inv(BHBt) operator");

   if(not allPassed || verbosity>=10) 
      os << ss.str(); 

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
