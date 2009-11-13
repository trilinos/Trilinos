#include "tLumping.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// PB-Package includes
#include "PB_Utilities.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

namespace PB {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

void tLumping::initializeTest()
{
   tolerance_ = 1.0e-15;
}

int tLumping::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLumping";

   status = test_lumping(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"lumping\" ... PASSED","   \"lumping\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_invLumping(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"invLumping\" ... PASSED","   \"invLumping\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tLumping...PASSED","tLumping...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tLumping...FAILED");
   }

   return failcount;
}

bool tLumping::test_lumping(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(tolerance_);

   Epetra_Map map(100,0,*GetComm());

   // A matrix...to be lumped
   Epetra_CrsMatrix A(Copy,map,5,false);
   int indices[5];
   double values[5] = {1,2,3,4,5};
   for(int i=0;i<A.NumMyRows()-5;i++) {
      int index = A.RowMap().GID(i);
      for(int j=0;j<5;j++)
         indices[j] = A.RowMap().GID(i+j);
      A.InsertGlobalValues(index,5,values,indices);
   }
   for(int i=A.NumMyRows()-5;i<A.NumMyRows();i++) {
      int index = A.RowMap().GID(i);
      for(int j=0;j<5;j++)
         indices[j] = A.RowMap().GID(j);
      A.InsertGlobalValues(index,5,values,indices);
   }
   A.FillComplete();

   // B matrix...already lumped
   Epetra_CrsMatrix B(Copy,map,1,true);
   double number = 15.0;
   for(int i=0;i<B.NumMyRows();i++) {
      int index = B.RowMap().GID(i);
      B.InsertGlobalValues(index,1,&number,&index);
   }
   B.FillComplete();

   PB::LinearOp pA = Thyra::epetraLinearOp(rcpFromRef(A));
   PB::LinearOp pB = Thyra::epetraLinearOp(rcpFromRef(B));
   PB::LinearOp lumpedA = getLumpedMatrix(pA);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *pB, *lumpedA, &fos );
      TEST_ASSERT(result,
             std::endl << "   tLummping::test_lumping "
             << ": Testing basic lumping functionality");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   return allPassed;
}

bool tLumping::test_invLumping(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(tolerance_);

   Epetra_Map map(100,0,*GetComm());

   // A matrix...to be lumped
   Epetra_CrsMatrix A(Copy,map,5,false);
   int indices[5];
   double values[5] = {1,2,3,4,5};
   for(int i=0;i<A.NumMyRows()-5;i++) {
      int index = A.RowMap().GID(i);
      for(int j=0;j<5;j++)
         indices[j] = A.RowMap().GID(i+j);
      A.InsertGlobalValues(index,5,values,indices);
   }
   for(int i=A.NumMyRows()-5;i<A.NumMyRows();i++) {
      int index = A.RowMap().GID(i);
      for(int j=0;j<5;j++)
         indices[j] = A.RowMap().GID(j);
      A.InsertGlobalValues(index,5,values,indices);
   }
   A.FillComplete();

   // B matrix...already lumped
   Epetra_CrsMatrix B(Copy,map,1,true);
   double number = 1.0/15.0;
   for(int i=0;i<B.NumMyRows();i++) {
      int index = B.RowMap().GID(i);
      B.InsertGlobalValues(index,1,&number,&index);
   }
   B.FillComplete();

   PB::LinearOp pA = Thyra::epetraLinearOp(rcpFromRef(A));
   PB::LinearOp pB = Thyra::epetraLinearOp(rcpFromRef(B));
   PB::LinearOp lumpedA = getInvLumpedMatrix(pA);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *pB, *lumpedA, &fos );
      TEST_ASSERT(result,
             std::endl << "   tLummping::test_invLumping "
             << ": Testing basic inv lumping functionality");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace PB
