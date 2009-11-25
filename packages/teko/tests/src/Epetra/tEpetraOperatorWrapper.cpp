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

#include "tEpetraOperatorWrapper.hpp"


namespace Teko {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;

void tEpetraOperatorWrapper::initializeTest() {}

int tEpetraOperatorWrapper::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tEpetraOperatorWrapper";

   status = test_functionality(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"functionality\" ... PASSED","   \"functionality\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tEpetraOperatorWrapper...PASSED","tEpetraOperatorWrapper...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tEpetraOperatorWrapper...FAILED");
   }

   return failcount;
}


bool tEpetraOperatorWrapper::test_functionality(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm = *GetComm();

   TEST_MSG("\n   tEpetraOperatorWrapper::test_functionality: "
         << "Running on " << comm.NumProc() << " processors");

   int nx = 39; // essentially random values
   int ny = 53;
   
   TEST_MSG( "   tEpetraOperatorWrapper::test_functionality: "
         << "Using Trilinos_Util to create test matrices" );

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> F = rcp(FGallery.GetMatrix(),false);

   Trilinos_Util::CrsMatrixGallery CGallery("laplace_2d",comm);
   CGallery.Set("nx",nx);
   CGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> C = rcp(CGallery.GetMatrix(),false);

   Trilinos_Util::CrsMatrixGallery BGallery("diag",comm);
   BGallery.Set("nx",nx*ny);
   BGallery.Set("a",5.0);
   RCP<Epetra_CrsMatrix> B = rcp(BGallery.GetMatrix(),false);

   Trilinos_Util::CrsMatrixGallery BtGallery("diag",comm);
   BtGallery.Set("nx",nx*ny);
   BtGallery.Set("a",3.0);
   RCP<Epetra_CrsMatrix> Bt = rcp(BtGallery.GetMatrix(),false);

   // load'em up in a thyra operator
   TEST_MSG("   tEpetraOperatorWrapper::test_functionality: "
         << " Building block2x2 Thyra matrix ... wrapping in EpetraOperatorWrapper");
   const RCP<const LinearOpBase<double> > A = Thyra::block2x2<double>(Thyra::epetraLinearOp(F),
                                                                      Thyra::epetraLinearOp(Bt),
                                                                      Thyra::epetraLinearOp(B),
                                                                      Thyra::epetraLinearOp(C),"A");

   const RCP<Thyra::EpetraOperatorWrapper> epetra_A = rcp(new Thyra::EpetraOperatorWrapper(A));

   // begin the tests!
   const Epetra_Map & rangeMap  = epetra_A->OperatorRangeMap();
   const Epetra_Map & domainMap = epetra_A->OperatorDomainMap();

   // check to see that the number of global elements is correct
   TEST_EQUALITY(rangeMap.NumGlobalElements(),2*nx*ny, 
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << "checking rangeMap size "
         << "( map = " << rangeMap.NumGlobalElements() 
         << ", true = " << 2*nx*ny << " )");

   // check to see that the number of global elements is correct
   TEST_EQUALITY(domainMap.NumGlobalElements(),2*nx*ny,
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << "checking domainMap size "
         << "( map = " << domainMap.NumGlobalElements() 
         << ", true = " << 2*nx*ny << " )" );

   // largest global ID should be one less then the # of elements
   TEST_EQUALITY(rangeMap.NumGlobalElements()-1,rangeMap.MaxAllGID(),
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking largest range element " 
         << "( largest = " << rangeMap.MaxAllGID() 
         << ", true = " <<  rangeMap.NumGlobalElements()-1 << " )" );
   TEST_EQUALITY(domainMap.NumGlobalElements()-1,domainMap.MaxAllGID(),
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking largest domain element " 
         << "( largest = " << domainMap.MaxAllGID() 
         << ", true = " <<  domainMap.NumGlobalElements()-1 << " )" );

   // create a vector to test: copyThyraIntoEpetra
   //////////////////////////////////////////////////////////////
   {
      const RCP<VectorBase<double> > tv = Thyra::createMember(A->domain());
      Thyra::randomize(-100.0,100.0,&*tv);
      const Thyra::ConstVector<double> handle_tv(tv);
      const RCP<const VectorBase<double> > tv_0 = Thyra::productVectorBase<double>(tv)->getVectorBlock(0);
      const RCP<const VectorBase<double> > tv_1 = Thyra::productVectorBase<double>(tv)->getVectorBlock(1);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_0(tv_0);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_1(tv_1);

      int off_0 = vv_0.globalOffset();
      int off_1 = vv_1.globalOffset();
      
      // create its Epetra counter part
      const RCP<Epetra_Vector> ev = rcp(new Epetra_Vector(epetra_A->OperatorDomainMap()));
      epetra_A->copyThyraIntoEpetra(handle_tv,*ev);

      // compare handle_tv to ev!
      TEST_EQUALITY(Thyra::dim(handle_tv),ev->GlobalLength(),
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking ThyraIntoEpetra copy "
         << "( thyra dim = " << Thyra::dim(handle_tv)
         << ", global length = " <<  ev->GlobalLength() << " )" );
      int numMyElements = domainMap.NumMyElements();
      bool compareThyraToEpetraValue = true;
      double tval = 0.0;
      for(int i=0;i<numMyElements;i++) {
         int gid = domainMap.GID(i);
         if(gid<nx*ny)
            tval = vv_0[gid-off_0];
         else
            tval = vv_1[gid-off_1-nx*ny];
         compareThyraToEpetraValue &= ((*ev)[i]==tval);
      }
      TEST_ASSERT(compareThyraToEpetraValue,
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " comparing Thyra to Epetra values" );
   }

   // create a vector to test: copyEpetraIntoThyra
   //////////////////////////////////////////////////////////////
   {
      // create an Epetra vector
      const RCP<Epetra_Vector> ev = rcp(new Epetra_Vector(epetra_A->OperatorDomainMap()));
      ev->Random();

      // create its thyra counterpart
      const RCP<VectorBase<double> > tv = Thyra::createMember(A->domain());
      const Thyra::Vector<double> handle_tv(tv);
      const RCP<const VectorBase<double> > tv_0 = Thyra::productVectorBase<double>(tv)->getVectorBlock(0);
      const RCP<const VectorBase<double> > tv_1 = Thyra::productVectorBase<double>(tv)->getVectorBlock(1);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_0(tv_0);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_1(tv_1);

      int off_0 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(tv_0->space())->localOffset();
      int off_1 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(tv_1->space())->localOffset();

      epetra_A->copyEpetraIntoThyra(*ev,handle_tv);
   
      // compare handle_tv to ev!
      TEST_EQUALITY(Thyra::dim(handle_tv),ev->GlobalLength(),
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking EpetraIntoThyra copy "
         << "( thyra dim = " << Thyra::dim(handle_tv)
         << ", global length = " <<  ev->GlobalLength() << " )" );
      int numMyElements = domainMap.NumMyElements();
      bool compareEpetraToThyraValue = true;
      double tval = 0.0;
      for(int i=0;i<numMyElements;i++) {
         int gid = domainMap.GID(i);
         if(gid<nx*ny)
            tval = vv_0[gid-off_0];
         else
            tval = vv_1[gid-off_1-nx*ny];
         compareEpetraToThyraValue &= ((*ev)[i]==tval);
      }
      TEST_ASSERT(compareEpetraToThyraValue,
            "   tEpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " comparing Thyra to Epetra values" );
   }

   return allPassed;
}


} // end Test namespace
} // end Teko namespace
