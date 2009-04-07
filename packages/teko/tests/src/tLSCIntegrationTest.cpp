#include "tLSCIntegrationTest.hpp"

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

// PB-Package includes
#include "Epetra/PB_EpetraHelpers.hpp"
#include "Epetra/PB_EpetraLSCHelpers.hpp"
#include "Epetra/PB_EpetraBlockPreconditioner.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"

#include "Thyra_EpetraLinearOp.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Test-rig
#include "Test_Utils.hpp"

#include "PB_InverseFactory.hpp"

namespace PB {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::epetraLinearOp;

void tLSCIntegrationTest::initializeTest()
{
   tolerance_ = 1.0e-7;

   velMap_  = rcp(new Epetra_Map(5890,0,*GetComm())); // map of velocity space
   prsMap_  = rcp(new Epetra_Map(769,0,*GetComm())); // map of pressure space
   fullMap_ = rcp(new Epetra_Map(769+5890,0,*GetComm())); // map of pressure space
}

void tLSCIntegrationTest::solveList(Teuchos::ParameterList & paramList,int vcycles)
{
   paramList.set("Linear Solver Type","AztecOO");
   paramList.sublist("Linear Solver Types")
            .sublist("AztecOO").sublist("Forward Solve").set("Max Iterations",500);
   paramList.sublist("Linear Solver Types")
            .sublist("AztecOO").sublist("Forward Solve").set("Tolerance",1e-12);
   paramList.sublist("Linear Solver Types")
            .sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","GMRES");
   paramList.sublist("Linear Solver Types")
            .sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Size of Krylov Subspace",500);
   paramList.sublist("Linear Solver Types")
            .sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",0);
   paramList.set("Preconditioner Type","ML");
   Teuchos::ParameterList & MLList = paramList.sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");

   // set default values for smoothed aggregation in MLList
   ML_Epetra::SetDefaults("SA",MLList);
   MLList.set("max levels",6);
   MLList.set("cycle applications",vcycles);
   MLList.set("coarse: type","Amesos-KLU");
}

void tLSCIntegrationTest::loadStableSystem()
{
   Epetra_CrsMatrix *F=0, *B=0, *Bt=0,*Qu=0;

   F=0; B=0; Bt=0; Qu=0;

   // read in stable discretization
   TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToCrsMatrix("./data/lsc_F_2.mm",*velMap_,*velMap_,*velMap_,F));
   TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToCrsMatrix("./data/lsc_B_2.mm",*prsMap_,*prsMap_,*velMap_,B));
   TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToCrsMatrix("./data/lsc_Bt_2.mm",*velMap_,*velMap_,*prsMap_,Bt));
   TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToCrsMatrix("./data/lsc_Qu_2.mm",*velMap_,*velMap_,*velMap_,Qu));

   sA_ = rcp(PB::Epetra::block2x2(F,Bt,B,0,"A"));

   // set stable matrix pointers
   sF_  = rcp(F); sB_  = rcp(B); sBt_ = rcp(Bt); sQu_ = rcp(Qu);

   // build an exporter to work around issue with MMFileToVector
   Epetra_Export exporter(*fullMap_,sA_->OperatorRangeMap());

   // read in RHS vector
   {
      Epetra_Vector *vfull=0, *temp=0;
 
      // read in rhs file 
      TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToVector("./data/lsc_rhs.mm",*fullMap_,vfull));

      // MMFileToVector is immplemented incompletely...thats why an exporter is used
      temp = new Epetra_Vector(sA_->OperatorRangeMap());
      temp->Export(*vfull,exporter,Insert);
      rhs_ = rcp(temp);

      delete vfull;
   }

   // read in solution vector
   {
      Epetra_Vector *vfull=0, *temp=0;
 
      // read in exact solution file 
      TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToVector("./data/lsc_exact_2.mm",*fullMap_,vfull));

      // MMFileToVector is immplemented incompletely...thats why an exporter is used
      temp = new Epetra_Vector(sA_->OperatorRangeMap());
      temp->Export(*vfull,exporter,Insert);
      sExact_ = rcp(temp);

      delete vfull;
   }
}

int tLSCIntegrationTest::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLSCIntegrationTest";

   status = test_withmassStable(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"withmassStable\" ... PASSED","   \"withmassStable\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_nomassStable(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"nomassStable\" ... PASSED","   \"nomassStable\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tLSCIntegrationTest...PASSED","tLSCIntegrationTest...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tLSCIntegrationTest...FAILED");
   }

   return failcount;
}

bool tLSCIntegrationTest::test_withmassStable(int verbosity,std::ostream & os)
{
   Teuchos::ParameterList paramList;
   solveList(paramList,8);

   RCP<const PB::InverseFactory> invFact = PB::invFactoryFromParamList(paramList,"ML");
   TEUCHOS_ASSERT(invFact!=Teuchos::null);

   bool status = false;
   bool allPassed = true;

   // load everything
   loadStableSystem();

   // if you get here you automatically pass the first test
   if(verbosity>=10 ) {
      os << std::endl << "   tLSCIntegrationTest::test_withmassStable: loading system ... " 
         << toString(true) << std::endl;
   }

   LinearOp Qu = epetraLinearOp(sQu_);
   const RCP<PB::NS::LSCStrategy> strategy = rcp(new PB::NS::InvLSCStrategy(invFact,Qu));
   const RCP<PB::BlockPreconditionerFactory> precFact = rcp(new PB::NS::LSCPreconditionerFactory(strategy));
   const RCP<PB::Epetra::EpetraBlockPreconditioner> prec = rcp(new PB::Epetra::EpetraBlockPreconditioner(precFact));
   prec->buildPreconditioner(*sA_);

   // B. Build solver and solve system
   Epetra_Vector x(sA_->OperatorDomainMap());
   x.Scale(0.0);

   // build Epetra problem
   Epetra_LinearProblem problem(&*sA_,&x,&*rhs_); // this doesn't take const arguments!

   AztecOO::AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres);
   solver.SetAztecOption(AZ_precond,AZ_none);
   solver.SetAztecOption(AZ_kspace,50);
   solver.SetAztecOption(AZ_output,AZ_none);
   solver.SetPrecOperator(&*prec);

   solver.Iterate(1000,1e-8);

   // check iteration count
   status = (solver.NumIters()<=16);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest::test_withmassStable " << toString(status) 
                      << ": # of iterations = " << solver.NumIters() << std::endl;
   }
   allPassed &= status;
 
   // check exact answer (versus IFISS solution)
   x.Update(-1.0,*sExact_,1.0); // x = x - x*
   double errnorm,exactnorm,relerr;
   x.Norm2(&errnorm);
   sExact_->Norm2(&exactnorm);
   status = ((relerr = errnorm/exactnorm) <= tolerance_);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest::test_withmassStable " << toString(status) 
                      << ": error in solution = " << std::scientific << relerr << " <= " << tolerance_ << std::endl;
   }
   allPassed &= status;

   return allPassed;
}

bool tLSCIntegrationTest::test_nomassStable(int verbosity,std::ostream & os)
{
   Teuchos::ParameterList paramList;
   solveList(paramList,8);

   RCP<const PB::InverseFactory> invFact = PB::invFactoryFromParamList(paramList,"ML");
   TEUCHOS_ASSERT(invFact!=Teuchos::null);

   bool status = false;
   bool allPassed = true;

   int vcycles = 8;

   // load everything
   loadStableSystem();

   // if you get here you automatically pass!
   if(verbosity>=10 ) {
      os << std::endl << "   tLSCIntegrationTest::test_nomassStable: loading system ... " 
         << toString(true) << std::endl;
   }

   const RCP<PB::NS::LSCStrategy> strategy = rcp(new PB::NS::InvLSCStrategy(invFact));
   const RCP<PB::BlockPreconditionerFactory> precFact = rcp(new PB::NS::LSCPreconditionerFactory(strategy));
   const RCP<PB::Epetra::EpetraBlockPreconditioner> prec = rcp(new PB::Epetra::EpetraBlockPreconditioner(precFact));
   prec->buildPreconditioner(*sA_);

   // B. Build solver and solve system
   Epetra_Vector x(sA_->OperatorDomainMap());
   x.Scale(0.0);

   // build Epetra problem
   Epetra_LinearProblem problem(&*sA_,&x,&*rhs_); // this doesn't take const arguments!

   AztecOO::AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres);
   solver.SetAztecOption(AZ_precond,AZ_none);
   solver.SetAztecOption(AZ_kspace,50);
   solver.SetAztecOption(AZ_output,AZ_none);
   solver.SetPrecOperator(&*prec);

   solver.Iterate(1000,1e-8);

   // check iteration count
   status = (solver.NumIters()==43);
   if(not status || verbosity>=10 ) { 
      os << std::endl << "   tLSCIntegrationTest::test_nomassStable " << toString(status) 
                      << ": # of iterations = " << solver.NumIters() << std::endl;
   }
   allPassed &= status;
 
   // check exact answer (versus IFISS solution)
   x.Update(-1.0,*sExact_,1.0); // x = x - x*
   double errnorm,exactnorm,relerr;
   x.Norm2(&errnorm);
   sExact_->Norm2(&exactnorm);
   status = ((relerr = errnorm/exactnorm) <= tolerance_);
   if(not status || verbosity==10 ) { 
      os << std::endl << "   tLSCIntegrationTest::test_nomassStable " << toString(status) 
                      << ": error in solution = " << std::scientific << relerr << std::endl;
   }
   allPassed &= status;

   return allPassed;
}

} // end namespace Tests
} // end namespace PB
