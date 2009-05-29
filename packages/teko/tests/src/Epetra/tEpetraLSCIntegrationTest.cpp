#include "tEpetraLSCIntegrationTest.hpp"

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

// Test-rig
#include "Test_Utils.hpp"

namespace PB {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::epetraLinearOp;

void tEpetraLSCIntegrationTest::initializeTest()
{
   tolerance_ = 1.0e-7;

   velMap_  = rcp(new Epetra_Map(5890,0,*GetComm())); // map of velocity space
   prsMap_  = rcp(new Epetra_Map(769,0,*GetComm())); // map of pressure space
   fullMap_ = rcp(new Epetra_Map(769+5890,0,*GetComm())); // map of pressure space
}

// construct a Epetra_Operator from an ML preconditioner
const RCP<const Epetra_Operator> tEpetraLSCIntegrationTest::mlSolve(const RCP<const Epetra_RowMatrix> mat,int cycles)
{
   // This is ripped directly from the ML user guide
   Teuchos::ParameterList MLList;

   // set default values for smoothed aggregation in MLList
   ML_Epetra::SetDefaults("SA",MLList);

   // overwrite with user’s defined parameters
   MLList.set("max levels",6);
   MLList.set("cycle applications",cycles);
   MLList.set("coarse: type","Amesos-KLU");

   // build preconditioner
   return rcp(new ML_Epetra::MultiLevelPreconditioner(*mat, MLList, true));
}

void tEpetraLSCIntegrationTest::loadStableSystem()
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

int tEpetraLSCIntegrationTest::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tEpetraLSCIntegrationTest";

   status = test_withmassStable(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"withmassStable\" ... PASSED","   \"withmassStable\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_nomassStable(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"nomassStable\" ... PASSED","   \"nomassStable\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tEpetraLSCIntegrationTest...PASSED","tEpetraLSCIntegrationTest...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tEpetraLSCIntegrationTest...FAILED");
   }

   return failcount;
}

bool tEpetraLSCIntegrationTest::test_withmassStable(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   int vcycles = 8;

   // load everything
   loadStableSystem();

   // if you get here you automatically pass!
   if(verbosity>=10 ) {
      os << std::endl << "   tEpetraLSCIntegrationTest::test_withmassStable: loading system ... " 
         << toString(true) << std::endl;
   }

   RCP<const Epetra_CrsMatrix> BBt;
   RCP<const Epetra_Vector> aiD; // should be null!

   PB::Epetra::buildLSCOperators(*sA_,*sQu_,BBt,aiD);

   TEUCHOS_ASSERT(aiD==Teuchos::null);

   const RCP<const Epetra_Operator> invF = mlSolve(sF_,vcycles);
   const RCP<const Epetra_Operator> invS = mlSolve(BBt,vcycles);

   // build inverse mass matrix
   const RCP<Epetra_Vector> invMass = rcp(new Epetra_Vector(*velMap_));
   TEST_FOR_EXCEPT(sQu_->ExtractDiagonalCopy(*invMass));
   TEST_FOR_EXCEPT(invMass->Reciprocal(*invMass));

   const RCP<PB::BlockPreconditionerFactory> precFact 
         = rcp(new PB::NS::LSCPreconditionerFactory(epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(&*invF))),
                                                    epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(&*invS))),
                                                    PB::Epetra::thyraDiagOp(invMass,invF->OperatorRangeMap())));
   const RCP<PB::Epetra::EpetraBlockPreconditioner> prec 
         = rcp(new PB::Epetra::EpetraBlockPreconditioner(precFact));
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
      os << std::endl << "   tEpetraLSCIntegrationTest::test_withmassStable " << toString(status) 
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
      os << std::endl << "   tEpetraLSCIntegrationTest::test_withmassStable " << toString(status) 
                      << ": error in solution = " << std::scientific << relerr << std::endl;
   }
   allPassed &= status;

   return allPassed;
}

bool tEpetraLSCIntegrationTest::test_nomassStable(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   int vcycles = 8;

   // load everything
   loadStableSystem();

   // if you get here you automatically pass!
   if(verbosity>=10 ) {
      os << std::endl << "   tEpetraLSCIntegrationTest::test_nomassStable: loading system ... " 
         << toString(true) << std::endl;
   }

   RCP<const Epetra_CrsMatrix> BBt;
   RCP<const Epetra_Vector> aiD; // should be null!

   PB::Epetra::buildLSCOperators(*sA_,BBt,aiD);

   TEUCHOS_ASSERT(aiD==Teuchos::null);

   const RCP<const Epetra_Operator> invF = mlSolve(sF_,vcycles);
   const RCP<const Epetra_Operator> invS = mlSolve(BBt,vcycles);

   const RCP<PB::BlockPreconditionerFactory> precFact 
         = rcp(new PB::NS::LSCPreconditionerFactory(epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(&*invF))),
                                                    epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(&*invS))),Teuchos::null));
   const RCP<PB::Epetra::EpetraBlockPreconditioner> prec 
         = rcp(new PB::Epetra::EpetraBlockPreconditioner(precFact));
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
      os << std::endl << "   tEpetraLSCIntegrationTest::test_nomassStable " << toString(status) 
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
      os << std::endl << "   tEpetraLSCIntegrationTest::test_nomassStable " << toString(status) 
                      << ": error in solution = " << std::scientific << relerr << std::endl;
   }
   allPassed &= status;

   return allPassed;
}

} // end namespace Tests
} // end namespace PB
