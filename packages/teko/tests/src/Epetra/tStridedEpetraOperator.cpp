// Thyra testing tools
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

// Thyra includes
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "tStridedEpetraOperator.hpp"

#include "Epetra/PB_StridedEpetraOperator.hpp"


namespace PB {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;

void tStridedEpetraOperator::initializeTest() 
{
   tolerance_ = 1e-10;
}

int tStridedEpetraOperator::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tStridedEpetraOperator";

   status = test_numvars_constr(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"numvars_constr\" ... PASSED","   \"numvars_constr\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_vector_constr(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"vector_constr\" ... PASSED","   \"vector_constr\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++; 

   status = test_reorder(verbosity,failstrm,0);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"reorder(flat reorder)\" ... PASSED","   \"reorder(flat reorder)\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_reorder(verbosity,failstrm,1);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"reorder(composite reorder = " << 1 
                      << ")\" ... PASSED","   \"reorder(composite reorder)\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_reorder(verbosity,failstrm,2);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"reorder(composite reorder = " << 2 
                      << ")\" ... PASSED","   \"reorder(composite reorder)\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tStridedEpetraOperator...PASSED","tStridedEpetraOperator...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tStridedEpetraOperator...FAILED");
   }

   return failcount;
}

bool tStridedEpetraOperator::test_numvars_constr(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm = *GetComm();

   TEST_MSG("\n   tStridedEpetraOperator::test_numvars: "
         << "Running on " << comm.NumProc() << " processors");

   // pick 
   int nx = 3 * 25 * comm.NumProc();
   int ny = 3 * 50 * comm.NumProc();


   // create a big matrix to play with
   // note: this matrix is not really strided
   //       however, I just need a nontrivial
   //       matrix to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> A = rcp(FGallery.GetMatrix(),false);

   int vars = 3;
   int width = 3;
   Epetra_MultiVector x(A->OperatorDomainMap(),width);
   Epetra_MultiVector ys(A->OperatorRangeMap(),width);
   Epetra_MultiVector y(A->OperatorRangeMap(),width);

   PB::Epetra::StridedEpetraOperator shell(vars,A);

   // test the operator against a lot of random vectors
   int numtests = 100;
   double max = 0.0;
   double min = 1.0;
   for(int i=0;i<numtests;i++) {
      double norm[width];
      x.Random();

      shell.Apply(x,y);
      A->Apply(x,ys);

      Epetra_MultiVector e(y);
      e.Update(-1.0,ys,1.0);
      e.Norm2(norm);

      for(int j=0;j<width;j++) {
         max = max>norm[j] ? max : norm[j];
         min = min<norm[j] ? min : norm[j];
      }
   }
   TEST_ASSERT(max>=min,
         "   tStridedEpetraOperator::test_numvars_constr: " << toString(status) << ": "
      << "sanity checked - " << max << " >= " << min);
   TEST_ASSERT(max<=tolerance_,
         "   tStridedEpetraOperator::test_numvars_constr: " << toString(status) << ": "
      << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << tolerance_ << " )");

   return allPassed;
}

bool tStridedEpetraOperator::test_vector_constr(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm = *GetComm();

   TEST_MSG("\n   tStridedEpetraOperator::test_vector_constr: "
         << "Running on " << comm.NumProc() << " processors");

   // pick 
   int nx = 3 * 25 * comm.NumProc();
   int ny = 3 * 50 * comm.NumProc();


   // create a big matrix to play with
   // note: this matrix is not really strided
   //       however, I just need a nontrivial
   //       matrix to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> A = rcp(FGallery.GetMatrix(),false);

   int width = 3;
   Epetra_MultiVector x(A->OperatorDomainMap(),width);
   Epetra_MultiVector ys(A->OperatorRangeMap(),width);
   Epetra_MultiVector y(A->OperatorRangeMap(),width);

   std::vector<int> vars;
   vars.push_back(2);
   vars.push_back(1);
   PB::Epetra::StridedEpetraOperator shell(vars,A);

   // test the operator against a lot of random vectors
   int numtests = 100;
   double max = 0.0;
   double min = 1.0;
   for(int i=0;i<numtests;i++) {
      double norm[width];
      x.Random();

      shell.Apply(x,y);
      A->Apply(x,ys);

      Epetra_MultiVector e(y);
      e.Update(-1.0,ys,1.0);
      e.Norm2(norm);

      for(int j=0;j<width;j++) {
         max = max>norm[j] ? max : norm[j];
         min = min<norm[j] ? min : norm[j];
      }
   }
   TEST_ASSERT(max>=min,
         "   tStridedEpetraOperator::test_vector_constr: " << toString(status) << ": "
      << "sanity checked - " << max << " >= " << min);
   TEST_ASSERT(max<=tolerance_,
         "   tStridedEpetraOperator::test_vector_constr: " << toString(status) << ": "
      << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << tolerance_ << " )");

   return allPassed;
}

bool tStridedEpetraOperator::test_reorder(int verbosity,std::ostream & os,int total)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm = *GetComm();

   std::string tstr = total ? "(composite reorder)" : "(flat reorder)";

   TEST_MSG("\n   tStridedEpetraOperator::test_reorder" << tstr << ": "
         << "Running on " << comm.NumProc() << " processors");

   // pick 
   int nx = 3 * 25 * comm.NumProc();
   int ny = 3 * 50 * comm.NumProc();

   // create a big matrix to play with
   // note: this matrix is not really strided
   //       however, I just need a nontrivial
   //       matrix to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> A = rcp(FGallery.GetMatrix(),false);

   int width = 3;
   Epetra_MultiVector x(A->OperatorDomainMap(),width);
   Epetra_MultiVector yf(A->OperatorRangeMap(),width);
   Epetra_MultiVector yr(A->OperatorRangeMap(),width);

   PB::Epetra::StridedEpetraOperator flatShell(3,A,"Af");
   PB::Epetra::StridedEpetraOperator reorderShell(3,A,"Ar");
 
   PB::BlockReorderManager brm;
   switch (total) {
   case 0:
      brm.SetNumBlocks(3);
      brm.SetBlock(0,1);
      brm.SetBlock(1,0);
      brm.SetBlock(2,2);
      break;
   case 1:
      brm.SetNumBlocks(2);
      brm.SetBlock(0,1);
      brm.GetBlock(1)->SetNumBlocks(2);
      brm.GetBlock(1)->SetBlock(0,0);
      brm.GetBlock(1)->SetBlock(1,2);
      break;
   case 2:
      brm.SetNumBlocks(2);
      brm.GetBlock(0)->SetNumBlocks(2);
      brm.GetBlock(0)->SetBlock(0,0);
      brm.GetBlock(0)->SetBlock(1,2);
      brm.SetBlock(1,1);
      break;
   }
   reorderShell.Reorder(brm);
   TEST_MSG("\n   tStridedEpetraOperator::test_reorder" << tstr << ": patern = " << brm.toString());

   TEST_MSG("\n   tStridedEpetraOperator::test_reorder" << tstr << ":\n");
   TEST_MSG("\n      " << Teuchos::describe(*reorderShell.getThyraOp(), Teuchos::VERB_HIGH)  << std::endl);

   // test the operator against a lot of random vectors
   int numtests = 10;
   double max = 0.0;
   double min = 1.0;
   for(int i=0;i<numtests;i++) {
      double norm[width];
      x.Random();

      flatShell.Apply(x,yf);
      reorderShell.Apply(x,yr);

      yf.Norm2(norm);
      yr.Norm2(norm);

      Epetra_MultiVector e(yf);
      e.Update(-1.0,yr,1.0);
      e.Norm2(norm);

      for(int j=0;j<width;j++) {
         max = max>norm[j] ? max : norm[j];
         min = min<norm[j] ? min : norm[j];
      }
   }
   TEST_ASSERT(max>=min,
         "   tStridedEpetraOperator::test_reorder" << tstr << ": " << toString(status) << ": "
      << "sanity checked - " << max << " >= " << min);
   TEST_ASSERT(max<=tolerance_,
         "   tStridedEpetraOperator::test_reorder" << tstr << ": "<< toString(status) << ": "
      << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << tolerance_ << " )");

   return allPassed;
}

} // end Test namespace
} // end PB namespace
