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

#include "tLSCIntegrationTest_tpetra.hpp"

#include <string>

// ML includes
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

// Teko-Package includes
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Test-rig
#include "Test_Utils.hpp"

// Tpetra includes
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include "Trilinos_Util_CrsMatrixGallery.h"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

void tLSCIntegrationTest_tpetra::initializeTest()
{
   tolerance_ = 1.0e-6;

   velMap_  = rcp(new Tpetra::Map<LO,GO,NT>(5890,0,GetComm_tpetra())); // map of velocity space
   prsMap_  = rcp(new Tpetra::Map<LO,GO,NT>(769,0,GetComm_tpetra())); // map of pressure space
   fullMap_ = rcp(new Tpetra::Map<LO,GO,NT>(769+5890,0,GetComm_tpetra())); // map of pressure space
}

void tLSCIntegrationTest_tpetra::solveList(Teuchos::ParameterList & paramList,int vcycles)
{
   paramList.set("Linear Solver Type","Belos");
   paramList.sublist("Linear Solver Types").sublist("Belos").set("Solver Type","Block GMRES");
   Teuchos::ParameterList & gmresList = paramList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES");
   gmresList.set( "Num Blocks", 100);                          // Maximum number of blocks in Krylov factorization
   gmresList.set( "Block Size", 1 );                          // Blocksize to be used by iterative solver
   gmresList.set( "Maximum Iterations", 1000 );                // Maximum number of iterations allowed
   gmresList.set( "Maximum Restarts", 15 );                   // Maximum number of restarts allowed
   gmresList.set( "Convergence Tolerance", 1.0e-9 );         // Relative convergence tolerance requested
   paramList.set("Preconditioner Type","Ifpack2");

   Teuchos::ParameterList & MLList = paramList.sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");

   // set default values for smoothed aggregation in MLList
   ML_Epetra::SetDefaults("SA",MLList);
   MLList.set("max levels",6);
   MLList.set("cycle applications",vcycles);
   MLList.set("coarse: type","Amesos-KLU");
}

void tLSCIntegrationTest_tpetra::loadStableSystem()
{
   // read in stable discretization
   RCP<const Tpetra::Map<LO,GO,NT> > junkMap = Teuchos::null;
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > nF_ = Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<ST,LO,GO,NT> >::readSparseFile("./data/lsc_F_2.mm",velMap_,junkMap,velMap_,velMap_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > nB_ = Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<ST,LO,GO,NT> >::readSparseFile("./data/lsc_B_2.mm",prsMap_,junkMap,velMap_,prsMap_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > nBt_ = Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<ST,LO,GO,NT> >::readSparseFile("./data/lsc_Bt_2.mm",velMap_,junkMap,prsMap_,velMap_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > nQu_ = Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<ST,LO,GO,NT> >::readSparseFile("./data/lsc_Qu_2.mm",velMap_,junkMap,velMap_,velMap_);

   sF_ = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST,LO,GO,NT> >(nF_);
   sB_ = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST,LO,GO,NT> >(nB_);
   sBt_ = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST,LO,GO,NT> >(nBt_);
   sQu_ = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST,LO,GO,NT> >(nQu_);

   Teko::LinearOp C;
   Teko::LinearOp tA_ = Thyra::block2x2<ST>(Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sF_->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sF_->getDomainMap()),sF_),Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sBt_->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sBt_->getDomainMap()),sBt_),Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sB_->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sB_->getDomainMap()),sB_),C,"A");
   sA_ = rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(tA_));

   // build an exporter to work around issue with MMFileToVector
   Tpetra::Export<LO,GO,NT> exporter(fullMap_,sA_->getRangeMap());

   // read in RHS vector
   {
      RCP<Tpetra::Vector<ST,LO,GO,NT> > vfull, temp;
 
      // read in rhs file 
      vfull = Tpetra::MatrixMarket::Reader<Tpetra::Vector<ST,LO,GO,NT> >::readVectorFile("./data/lsc_rhs.mm",GetComm_tpetra(),fullMap_);

      temp = rcp(new Tpetra::Vector<ST,LO,GO,NT>(sA_->getRangeMap()));
      temp->doExport(*vfull,exporter,Tpetra::INSERT);
      rhs_ = temp;
   }

   // read in solution vector
   {
      RCP<Tpetra::Vector<ST,LO,GO,NT> > vfull, temp;
 
      // read in exact solution file 
      vfull = Tpetra::MatrixMarket::Reader<Tpetra::Vector<ST,LO,GO,NT> >::readVectorFile("./data/lsc_exact_2.mm",GetComm_tpetra(),fullMap_);

      temp = rcp(new Tpetra::Vector<ST,LO,GO,NT>(sA_->getRangeMap()));
      temp->doExport(*vfull,exporter,Tpetra::INSERT);
      sExact_ = temp;
   }
}

int tLSCIntegrationTest_tpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLSCIntegrationTest_tpetra";

   status = test_withmassStable(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"withmassStable\" ... PASSED","   \"withmassStable\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_nomassStable(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"nomassStable\" ... PASSED","   \"nomassStable\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_plConstruction(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"plConstruction\" ... PASSED","   \"plConstruction\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tLSCIntegrationTest_tpetra...PASSED","tLSCIntegrationTest_tpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tLSCIntegrationTest_tpetra...FAILED");
   }

   return failcount;
}

bool tLSCIntegrationTest_tpetra::test_withmassStable(int verbosity,std::ostream & os)
{
   typedef Tpetra::MultiVector<ST,LO,GO,NT> MV;
   typedef Tpetra::Operator<ST,LO,GO,NT> OP;

   Teuchos::ParameterList paramList;
   solveList(paramList,8);

   RCP<Teko::InverseFactory> invFact = Teko::invFactoryFromParamList(paramList,"Belos");
   TEUCHOS_ASSERT(invFact!=Teuchos::null);

   bool status = false;
   bool allPassed = true;

   // load everything
   loadStableSystem();

   // if you get here you automatically pass the first test
   if(verbosity>=10 ) {
      os << std::endl << "   tLSCIntegrationTest_tpetra::test_withmassStable: loading system ... " 
         << toString(true) << std::endl;
   }

   LinearOp Qu = Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sQu_->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(sQu_->getDomainMap()),sQu_);
   const RCP<Teko::NS::LSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(invFact,Qu));
   const RCP<Teko::BlockPreconditionerFactory> precFact = rcp(new Teko::NS::LSCPreconditionerFactory(strategy));
   const RCP<Teko::TpetraHelpers::TpetraBlockPreconditioner> prec = rcp(new Teko::TpetraHelpers::TpetraBlockPreconditioner(precFact));
   prec->buildPreconditioner(sA_);

   // Build solver and solve system
   Tpetra::Vector<ST,LO,GO,NT> x(sA_->getDomainMap());
   x.scale(0.0);
   RCP<MV> x_mv = Teuchos::rcp_dynamic_cast<MV>(rcpFromRef(x));
   RCP<const MV> rhs_mv = Teuchos::rcp_dynamic_cast<const MV>(rhs_);

   // build Belos problem
   Belos::LinearProblem<ST,MV,OP> problem(sA_,x_mv,rhs_mv);
   problem.setRightPrec(prec);
   bool set = problem.setProblem();
   TEUCHOS_ASSERT(set);

   // build Belos solver
   RCP< Belos::SolverManager<ST,MV,OP> > solver
     = rcp( new Belos::BlockGmresSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&paramList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES"),false)));

   // do Belos solve
   solver->solve();

   // check iteration count
   int numIters = solver->getNumIters();
   status = (numIters<=18);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest_tpetra::test_withmassStable " << toString(status) 
                      << ": # of iterations = " << numIters << " (should be 18)" << std::endl;
   }
   allPassed &= status;
 
   // check exact answer (versus IFISS solution)
   x.update(-1.0,*sExact_,1.0); // x = x - x*
   ST errnorm,exactnorm,relerr;
   errnorm = x.norm2();
   exactnorm = sExact_->norm2();
   status = ((relerr = errnorm/exactnorm) <= tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest_tpetra::test_withmassStable " << toString(status) 
                      << ": error in solution = " << std::scientific << relerr << " <= " << tolerance_ << std::endl;
   }
   allPassed &= status;

   return allPassed;
}

bool tLSCIntegrationTest_tpetra::test_nomassStable(int verbosity,std::ostream & os)
{
   typedef Tpetra::MultiVector<ST,LO,GO,NT> MV;
   typedef Tpetra::Operator<ST,LO,GO,NT> OP;

   Teuchos::ParameterList paramList;
   solveList(paramList,8);

   RCP<Teko::InverseFactory> invFact = Teko::invFactoryFromParamList(paramList,"Belos");
   TEUCHOS_ASSERT(invFact!=Teuchos::null);

   bool status = false;
   bool allPassed = true;

   // int vcycles = 8;

   // load everything
   loadStableSystem();

   // if you get here you automatically pass!
   if(verbosity>=10 ) {
      os << std::endl << "   tLSCIntegrationTest_tpetra::test_nomassStable: loading system ... " 
         << toString(true) << std::endl;
   }

   const RCP<Teko::NS::LSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(invFact));
   const RCP<Teko::BlockPreconditionerFactory> precFact = rcp(new Teko::NS::LSCPreconditionerFactory(strategy));
   const RCP<Teko::TpetraHelpers::TpetraBlockPreconditioner> prec = rcp(new Teko::TpetraHelpers::TpetraBlockPreconditioner(precFact));
   prec->buildPreconditioner(sA_);

   // Build solver and solve system
   Tpetra::Vector<ST,LO,GO,NT> x(sA_->getDomainMap());
   x.scale(0.0);
   RCP<MV> x_mv = Teuchos::rcp_dynamic_cast<MV>(rcpFromRef(x));
   RCP<const MV> rhs_mv = Teuchos::rcp_dynamic_cast<const MV>(rhs_);

   // build Belos problem
   Belos::LinearProblem<ST,MV,OP> problem(sA_,x_mv,rhs_mv);
   problem.setRightPrec(prec);
   bool set = problem.setProblem();
   TEUCHOS_ASSERT(set);

   // build Belos solver
   RCP< Belos::SolverManager<ST,MV,OP> > solver
     = rcp( new Belos::BlockGmresSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&paramList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES"),false)));

   // do Belos solve
   solver->solve();

   // check iteration count
   int numIters = solver->getNumIters();
   status = (numIters<=30);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest_tpetra::test_nomassStable " << toString(status) 
                      << ": # of iterations = " << numIters << " (should be 30)" << std::endl;
   }
   allPassed &= status;
 
   // check exact answer (versus IFISS solution)
   x.update(-1.0,*sExact_,1.0); // x = x - x*
   ST errnorm,exactnorm,relerr;
   errnorm = x.norm2();
   exactnorm = sExact_->norm2();
   status = ((relerr = errnorm/exactnorm) <= tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest_tpetra::test_nomassStable " << toString(status) 
                      << ": error in solution = " << std::scientific << relerr << " <= " << tolerance_ << std::endl;
   }
   allPassed &= status;

   return allPassed;
}

bool tLSCIntegrationTest_tpetra::test_plConstruction(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   using Teuchos::ParameterList;

   RCP<Teko::PreconditionerFactory> precFact;
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

   /////////////////////////////////////////////////////////////////////////////

   ParameterList pl; 
   pl.set("Inverse Type", "Amesos");
   pl.set("Inverse Velocity Type", "Ifpack");
   pl.set("Inverse Pressure Type", "Ifpack");
   pl.set("Ignore Boundary Rows",true);
   pl.set("Use LDU",true);

   precFact = Teko::PreconditionerFactory::buildPreconditionerFactory("NS LSC", pl, invLib);

   TEST_ASSERT(precFact!=Teuchos::null,
         std::endl << "   tLSCIntegrationTest_tpetra::test_plConstruction " << toString(status)
                   << ": building \"Basic Inverse\" with out sublist");

   /////////////////////////////////////////////////////////////////////////////

   ParameterList parPl;
   parPl.set("Strategy Name","Basic Inverse");
   parPl.set("Strategy Settings",pl);

   precFact = Teko::PreconditionerFactory::buildPreconditionerFactory("NS LSC", parPl, invLib);

   TEST_ASSERT(precFact!=Teuchos::null,
         std::endl << "   tLSCIntegrationTest_tpetra::test_plConstruction " << toString(status)
                   << ": building \"Basic Inverse\" with sublist");

   /////////////////////////////////////////////////////////////////////////////

   try {
      parPl.set("Strategy Name","The Cat");
      precFact = Teko::PreconditionerFactory::buildPreconditionerFactory("NS LSC", parPl, invLib);

      TEST_ASSERT(false,
         std::endl << "   tLSCIntegrationTest_tpetra::test_plConstruction " << toString(status)
                   << ": using failing strategy to build LSC factory...exception expected");
   }
   catch (const std::runtime_error & e) { }

   /////////////////////////////////////////////////////////////////////////////

   try {
      pl.set("Strategy Name","The Cat");
      precFact = Teko::PreconditionerFactory::buildPreconditionerFactory("NS LSC", pl, invLib);

      TEST_ASSERT(false,
            std::endl << "   tLSCIntegrationTest_tpetra::test_plConstruction " << toString(status)
                      << ": using failing strategy to build LSC factory...exception expected");
   }
   catch (const std::runtime_error & e) { }

   /////////////////////////////////////////////////////////////////////////////

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
