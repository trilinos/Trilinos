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

// Thyra testing tools
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

// Thyra includes
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraHelpers.hpp"

#include "tTpetraOperatorWrapper.hpp"


namespace Teko {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;

void tTpetraOperatorWrapper::initializeTest() {}

int tTpetraOperatorWrapper::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tTpetraOperatorWrapper";

   status = test_functionality(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"functionality\" ... PASSED","   \"functionality\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tTpetraOperatorWrapper...PASSED","tTpetraOperatorWrapper...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tTpetraOperatorWrapper...FAILED");
   }

   return failcount;
}


bool tTpetraOperatorWrapper::test_functionality(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm_epetra = *GetComm();

   TEST_MSG("\n   tTpetraOperatorWrapper::test_functionality: "
         << "Running on " << comm_epetra.NumProc() << " processors");

   int nx = 39; // essentially random values
   int ny = 53;
   
   TEST_MSG( "   tTpetraOperatorWrapper::test_functionality: "
         << "Using Trilinos_Util to create test matrices" );

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm_epetra,false);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   Epetra_CrsMatrix & epetraF = FGallery.GetMatrixRef();
   RCP<const Tpetra::CrsMatrix<ST,LO,GO,NT> > tpetraF = Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraF),comm_tpetra);

   Trilinos_Util::CrsMatrixGallery CGallery("laplace_2d",comm_epetra,false);
   CGallery.Set("nx",nx);
   CGallery.Set("ny",ny);
   Epetra_CrsMatrix & epetraC = CGallery.GetMatrixRef();
   RCP<const Tpetra::CrsMatrix<ST,LO,GO,NT> > tpetraC = Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraC),comm_tpetra);

   Trilinos_Util::CrsMatrixGallery BGallery("diag",comm_epetra,false);
   BGallery.Set("nx",nx*ny);
   BGallery.Set("a",5.0);
   Epetra_CrsMatrix & epetraB = BGallery.GetMatrixRef();
   RCP<const Tpetra::CrsMatrix<ST,LO,GO,NT> > tpetraB = Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraB),comm_tpetra);

   Trilinos_Util::CrsMatrixGallery BtGallery("diag",comm_epetra,false);
   BtGallery.Set("nx",nx*ny);
   BtGallery.Set("a",3.0);
   Epetra_CrsMatrix & epetraBt = BtGallery.GetMatrixRef();
   RCP<const Tpetra::CrsMatrix<ST,LO,GO,NT> > tpetraBt = Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraBt),comm_tpetra);

   // load'em up in a thyra operator
   TEST_MSG("   tTpetraOperatorWrapper::test_functionality: "
         << " Building block2x2 Thyra matrix ... wrapping in TpetraOperatorWrapper");
   const RCP<const LinearOpBase<double> > A = Thyra::block2x2<double>(Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraF->getDomainMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraF->getRangeMap()),tpetraF),
                                                                      Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraBt->getDomainMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraBt->getRangeMap()),tpetraBt),
                                                                      Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraB->getDomainMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraB->getRangeMap()),tpetraB),
                                                                      Thyra::constTpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraC->getDomainMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraC->getRangeMap()),tpetraC),"A");

   // const RCP<Thyra::TpetraOperatorWrapper> epetra_A = rcp(new Thyra::TpetraOperatorWrapper(A));
   const RCP<Teko::TpetraHelpers::TpetraOperatorWrapper> tpetra_A = rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(A));

   // begin the tests!
   const RCP<const Tpetra::Map<LO,GO,NT> > & rangeMap  = tpetra_A->getRangeMap();
   const RCP<const Tpetra::Map<LO,GO,NT> >& domainMap = tpetra_A->getDomainMap();

   // check to see that the number of global elements is correct
   TEST_EQUALITY(rangeMap->getGlobalNumElements(),(Tpetra::global_size_t) 2*nx*ny, 
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << "checking rangeMap size "
         << "( map = " << rangeMap->getGlobalNumElements() 
         << ", true = " << 2*nx*ny << " )");

   // check to see that the number of global elements is correct
   TEST_EQUALITY(domainMap->getGlobalNumElements(),(Tpetra::global_size_t) 2*nx*ny,
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << "checking domainMap size "
         << "( map = " << domainMap->getGlobalNumElements() 
         << ", true = " << 2*nx*ny << " )" );

   // largest global ID should be one less then the # of elements
   TEST_EQUALITY(rangeMap->getGlobalNumElements()-1,(Tpetra::global_size_t) rangeMap->getMaxAllGlobalIndex(),
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking largest range element " 
         << "( largest = " << rangeMap->getMaxAllGlobalIndex() 
         << ", true = " <<  rangeMap->getGlobalNumElements()-1 << " )" );
   TEST_EQUALITY(domainMap->getGlobalNumElements()-1,(Tpetra::global_size_t) domainMap->getMaxAllGlobalIndex(),
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking largest domain element " 
         << "( largest = " << domainMap->getMaxAllGlobalIndex() 
         << ", true = " <<  domainMap->getGlobalNumElements()-1 << " )" );

   RCP<const Teko::TpetraHelpers::MappingStrategy> ms = tpetra_A->getMapStrategy();

   // create a vector to test: copyThyraIntoTpetra
   //////////////////////////////////////////////////////////////
   {
      const RCP<MultiVectorBase<ST> > tv = Thyra::createMembers(A->domain(),1);
      Thyra::randomize(-100.0,100.0,tv.ptr());
      // const Thyra::ConstVector<double> handle_tv(tv);
      const RCP<const MultiVectorBase<ST> > tv_0 
            = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> >(tv)->getMultiVectorBlock(0);
      const RCP<const MultiVectorBase<ST> > tv_1 
            = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> >(tv)->getMultiVectorBlock(1);
      const Thyra::ConstDetachedSpmdVectorView<ST> vv_0(tv_0->col(0));
      const Thyra::ConstDetachedSpmdVectorView<ST> vv_1(tv_1->col(0));

      LO off_0 = vv_0.globalOffset();
      LO off_1 = vv_1.globalOffset();
      
      // create its Tpetra counter part
      const RCP<Tpetra::Vector<ST,LO,GO,NT> > ev = rcp(new Tpetra::Vector<ST,LO,GO,NT>(tpetra_A->getDomainMap()));
      ms->copyThyraIntoTpetra(tv,*ev);

      // compare tv to ev!
      TEST_EQUALITY((Tpetra::global_size_t) tv->range()->dim(),ev->getGlobalLength(),
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking ThyraIntoTpetra copy "
         << "( thyra dim = " << tv->range()->dim()
         << ", global length = " <<  ev->getGlobalLength() << " )" );
      LO numMyElements = domainMap->getNodeNumElements();
TEST_MSG("domainMap->getNodeNumElements() = " << domainMap->getNodeNumElements());
      bool compareThyraToTpetraValue = true;
      ST tval = 0.0;
      for(LO i=0;i<numMyElements;i++) {
         GO gid = domainMap->getGlobalElement(i);
         if(gid-off_0 < nx*ny){
            tval = vv_0[gid-off_0];
         }else{
            tval = vv_1[gid-off_1-nx*ny];
         }
         compareThyraToTpetraValue &= (ev->get1dView()[i]==tval);
      }
      TEST_ASSERT(compareThyraToTpetraValue,
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " comparing Thyra to Tpetra values" );
   }

   // create a vector to test: copyTpetraIntoThyra
   //////////////////////////////////////////////////////////////
   {
      // create an Tpetra vector
      const RCP<Tpetra::Vector<ST,LO,GO,NT> > ev = rcp(new Tpetra::Vector<ST,LO,GO,NT>(tpetra_A->getDomainMap()));
      ev->randomize();

      // create its thyra counterpart
      const RCP<MultiVectorBase<ST> > tv = Thyra::createMembers(A->domain(),1);
      const RCP<const MultiVectorBase<ST> > tv_0 
            = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> >(tv)->getMultiVectorBlock(0);
      const RCP<const MultiVectorBase<ST> > tv_1 
            = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> >(tv)->getMultiVectorBlock(1);
      const Thyra::ConstDetachedSpmdVectorView<ST> vv_0(tv_0->col(0));
      const Thyra::ConstDetachedSpmdVectorView<ST> vv_1(tv_1->col(0));

      LO off_0 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<ST> >(tv_0->range())->localOffset();
      LO off_1 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<ST> >(tv_1->range())->localOffset();

      ms->copyTpetraIntoThyra(*ev,tv.ptr());
   
      // compare handle_tv to ev!
      TEST_EQUALITY((Tpetra::global_size_t) tv->range()->dim(),ev->getGlobalLength(),
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " checking TpetraIntoThyra copy "
         << "( thyra dim = " << tv->range()->dim()
         << ", global length = " <<  ev->getGlobalLength() << " )" );
      LO numMyElements = domainMap->getNodeNumElements();
      bool compareTpetraToThyraValue = true;
      ST tval = 0.0;
      for(LO i=0;i<numMyElements;i++) {
         GO gid = domainMap->getGlobalElement(i);
         if(gid-off_0 < nx*ny){
            tval = vv_0[gid-off_0];
         }else{
            tval = vv_1[gid-off_1-nx*ny];
         }
         compareTpetraToThyraValue &= (ev->get1dView()[i]==tval);
      }
      TEST_ASSERT(compareTpetraToThyraValue,
            "   tTpetraOperatorWrapper::test_functionality: " << toString(status) << ": "
         << " comparing Thyra to Tpetra values" );
   }

   return allPassed;
}


} // end Test namespace
} // end Teko namespace
