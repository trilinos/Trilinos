#include "tExplicitOps.hpp"

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

void tExplicitOps::initializeTest()
{
   const Epetra_Comm & comm = *GetComm();

   tolerance_ = 1.0e-4;

   int nx = 39; // essentially random values
   int ny = 53;

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   Epetra_CrsMatrix & epetraF = FGallery.GetMatrixRef();
   F_ = Thyra::nonconstEpetraLinearOp(rcp(new Epetra_CrsMatrix(epetraF)));

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery GGallery("laplace_2d",comm);
   GGallery.Set("nx",nx);
   GGallery.Set("ny",ny);
   Epetra_CrsMatrix & epetraG = GGallery.GetMatrixRef();
   G_ = Thyra::nonconstEpetraLinearOp(rcp(new Epetra_CrsMatrix(epetraG)));

   RCP<Epetra_Vector> v = rcp(new Epetra_Vector(epetraF.OperatorRangeMap()));
   v->Random();
   RCP<Thyra::VectorBase<double> > tV = Thyra::create_Vector(v,Thyra::create_VectorSpace(rcpFromRef(epetraF.RowMap()))); 
   D_ = Thyra::diagonal(tV);
}

int tExplicitOps::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tExplicitOps";

   status = test_mult_diagScaleMatProd(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"mult_diagScaleMatProd\" ... PASSED","   \"mult_diagScaleMatProd\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_mult_diagScaling(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"mult_diagScaling\" ... PASSED","   \"mult_diagScaling\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_mult_modScaleMatProd(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"mult_modScaleMatProd\" ... PASSED","   \"mult_modScaleMatProd\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tExplicitOps...PASSED","tExplicitOps...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tExplicitOps...FAILED");
   }

   return failcount;
}

bool tExplicitOps::test_mult_diagScaleMatProd(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<double> > thyOp;
   PB::LinearOp expOp;

   thyOp = Thyra::multiply(PB::scale(-4.0,F_),D_,PB::adjoint(G_));
   expOp = PB::explicitMultiply(PB::scale(-4.0,F_),D_,PB::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, &fos );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_diagScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   thyOp = PB::multiply(PB::scale(-4.0,F_),PB::adjoint(G_));
   expOp = PB::explicitMultiply(PB::scale(-4.0,F_),PB::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, &fos );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_diagScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps::test_mult_diagScaling(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<double> > thyOp;
   PB::LinearOp expOp;

   thyOp = PB::multiply(PB::scale(-4.0,F_),D_);
   expOp = PB::explicitMultiply(PB::scale(-4.0,F_),D_);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, &fos );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_diagScaleMatProd "
             << ": Testing diagonal scaling");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   thyOp = PB::multiply(D_,PB::scale(-9.0,F_));
   expOp = PB::explicitMultiply(D_,PB::scale(-9.0,F_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, &fos );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_diagScaleMatProd "
             << ": Testing diagonal scaling");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps::test_mult_modScaleMatProd(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);

   PB::LinearOp thyOp;
   PB::ModifiableLinearOp expOp;

   thyOp = PB::multiply(PB::scale(-4.0,F_),D_,PB::adjoint(G_));
   expOp = PB::explicitMultiply(PB::scale(-4.0,F_),D_,PB::adjoint(G_),expOp);

   RCP<const Epetra_Operator> eop1 = Thyra::get_Epetra_Operator(*expOp);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, &fos );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_modScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*F_))->Scale(5.0);
   Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*G_))->Scale(2.0);

   int numEntries = 0;
   double * values;

   // do some random violence (oh my brothers) to one row
   Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*F_))->ExtractMyRowView(3,numEntries,values);
   for(int i=0;i<numEntries;i++) values[i] *= values[i]*double(i+1)*0.92;

   // do some random violence (oh my brothers) to one row
   Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*F_))->ExtractMyRowView(7,numEntries,values);
   for(int i=0;i<numEntries;i++) values[i] *= values[i]*double(i+1)*0.92;

   // perform the next test
   thyOp = PB::multiply(PB::scale(-4.0,F_),D_,PB::adjoint(G_));
   expOp = PB::explicitMultiply(PB::scale(-4.0,F_),D_,PB::adjoint(G_),expOp);

   RCP<const Epetra_Operator> eop2 = Thyra::get_Epetra_Operator(*expOp);
   TEUCHOS_ASSERT(eop1==eop2);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, &fos );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_modScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace PB
