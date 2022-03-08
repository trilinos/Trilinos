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
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
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

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "tBlockedTpetraOperator.hpp"

#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraHelpers.hpp"

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

void tBlockedTpetraOperator::buildBlockGIDs(std::vector<std::vector<GO> > & gids,
                                            const Tpetra::Map<LO,GO,NT> & map) const
{
   LO numLocal = map.getLocalNumElements();
   LO numHalf = numLocal/2;
   numHalf += ((numHalf % 2 == 0) ? 0 : 1); 

   gids.clear();
   gids.resize(3);

   std::vector<GO> & blk0 = gids[0];
   std::vector<GO> & blk1 = gids[1];
   std::vector<GO> & blk2 = gids[2];
   
   // loop over global IDs: treat first block as strided
   GO gid = -1;
   for(LO i=0;i<numHalf;i+=2) {
      gid = map.getGlobalElement(i);
      blk0.push_back(gid);

      gid = map.getGlobalElement(i+1);
      blk1.push_back(gid);
   }

   // loop over global IDs: treat remainder as contiguous
   for(LO i=numHalf;i<numLocal;i++) {
      gid = map.getGlobalElement(i);
      blk2.push_back(gid);
   }

   TEUCHOS_ASSERT(LO(blk0.size()+blk1.size()+blk2.size())==numLocal);
}

void tBlockedTpetraOperator::initializeTest() 
{
   tolerance_ = 1e-14;
}

int tBlockedTpetraOperator::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tBlockedTpetraOperator";

   status = test_vector_constr(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"vector_constr\" ... PASSED","   \"vector_constr\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++; 

   status = test_reorder(verbosity,failstrm,0);
   Teko_TEST_MSG(stdstrm,1,"   \"reorder(flat reorder)\" ... PASSED","   \"reorder(flat reorder)\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_reorder(verbosity,failstrm,1);
   Teko_TEST_MSG(stdstrm,1,"   \"reorder(composite reorder = " << 1 
                      << ")\" ... PASSED","   \"reorder(composite reorder)\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_reorder(verbosity,failstrm,2);
   Teko_TEST_MSG(stdstrm,1,"   \"reorder(composite reorder = " << 2 
                      << ")\" ... PASSED","   \"reorder(composite reorder)\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tBlockedTpetraOperator...PASSED","tBlockedTpetraOperator...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tBlockedTpetraOperator...FAILED");
   }

   return failcount;
}

bool tBlockedTpetraOperator::test_vector_constr(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm_epetra = *GetComm();

   TEST_MSG("\n   tBlockedTpetraOperator::test_vector_constr: "
         << "Running on " << comm_epetra.NumProc() << " processors");

   // pick 
   int nx = 5;// * comm_epetra.NumProc();
   int ny = 5;// * comm_epetra.NumProc();

   // create a big matrix to play with
   // note: this matrix is not really strided
   //       however, I just need a nontrivial
   //       matrix to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm_epetra,false);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> epetraA = rcp(FGallery.GetMatrix(),false);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > A = Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(epetraA,comm_tpetra);
   ST beforeNorm = A->getFrobeniusNorm();

   int width = 3;
   Tpetra::MultiVector<ST,LO,GO,NT> x(A->getDomainMap(),width);
   Tpetra::MultiVector<ST,LO,GO,NT> ys(A->getRangeMap(),width);
   Tpetra::MultiVector<ST,LO,GO,NT> y(A->getRangeMap(),width);

   std::vector<std::vector<GO> > vars;
   buildBlockGIDs(vars,*A->getRowMap());

   Teko::TpetraHelpers::BlockedTpetraOperator shell(vars,A);

   // test the operator against a lot of random vectors
   int numtests = 50;
   ST max = 0.0;
   ST min = 1.0;
   for(int i=0;i<numtests;i++) {
      std::vector<ST> norm(width);
      std::vector<ST> rel(width);
      x.randomize();

      shell.apply(x,y);
      A->apply(x,ys);

      Tpetra::MultiVector<ST,LO,GO,NT> e(y,Teuchos::Copy);
      e.update(-1.0,ys,1.0);
      e.norm2(Teuchos::ArrayView<ST>(norm));

      // compute relative error
      ys.norm2(Teuchos::ArrayView<ST>(rel));
      for(int j=0;j<width;j++) {
         max = max>norm[j]/rel[j] ? max : norm[j]/rel[j];
         min = min<norm[j]/rel[j] ? min : norm[j]/rel[j];
      }
   }
   TEST_ASSERT(max>=min,
         "\n   tBlockedTpetraOperator::test_vector_constr: " << toString(status) << ": "
      << "sanity checked - " << max << " >= " << min);
   TEST_ASSERT(max<=tolerance_,
         "\n   tBlockedTpetraOperator::test_vector_constr: " << toString(status) << ": "     << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << tolerance_ << " )");

   A->scale(2.0); //double everything

   ST afterNorm = A->getFrobeniusNorm();
   TEST_ASSERT(beforeNorm!=afterNorm,
         "\n   tBlockedTpetraOperator::test_vector_constr " << toString(status) << ": "
      << "verify matrix has been modified");

   shell.RebuildOps();

   // test the operator against a lot of random vectors
   numtests = 50;
   max = 0.0;
   min = 1.0;
   for(int i=0;i<numtests;i++) {
      std::vector<ST> norm(width);
      std::vector<ST> rel(width);
      x.randomize();

      shell.apply(x,y);
      A->apply(x,ys);

      Tpetra::MultiVector<ST,LO,GO,NT> e(y,Teuchos::Copy);
      e.update(-1.0,ys,1.0);
      e.norm2(Teuchos::ArrayView<ST>(norm));

      // compute relative error
      ys.norm2(Teuchos::ArrayView<ST>(rel));
      for(int j=0;j<width;j++) {
         max = max>norm[j]/rel[j] ? max : norm[j]/rel[j];
         min = min<norm[j]/rel[j] ? min : norm[j]/rel[j];
      }
   }
   TEST_ASSERT(max>=min,
         "\n   tBlockedTpetraOperator::test_vector_constr (rebuild): " << toString(status) << ": "
      << "sanity checked - " << max << " >= " << min);
   TEST_ASSERT(max<=tolerance_,
         "\n   tBlockedTpetraOperator::test_vector_constr (rebuild): " << toString(status) << ": "
      << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << tolerance_ << " )");

   return allPassed;
}

bool tBlockedTpetraOperator::test_reorder(int verbosity,std::ostream & os,int total)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm_epetra = *GetComm();

   std::string tstr = total ? "(composite reorder)" : "(flat reorder)";

   TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ": "
         << "Running on " << comm_epetra.NumProc() << " processors");

   // pick 
   int nx = 5;//3 * 25 * comm_epetra.NumProc();
   int ny = 5;//3 * 50 * comm_epetra.NumProc();

   // create a big matrix to play with
   // note: this matrix is not really strided
   //       however, I just need a nontrivial
   //       matrix to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm_epetra,false);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> epetraA = rcp(FGallery.GetMatrix(),false);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > A = Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(epetraA,comm_tpetra);

   int width = 3;
   Tpetra::MultiVector<ST,LO,GO,NT> x(A->getDomainMap(),width);
   Tpetra::MultiVector<ST,LO,GO,NT> yf(A->getRangeMap(),width);
   Tpetra::MultiVector<ST,LO,GO,NT> yr(A->getRangeMap(),width);

   std::vector<std::vector<GO> > vars;
   buildBlockGIDs(vars,*A->getRowMap());

   Teko::TpetraHelpers::BlockedTpetraOperator flatShell(vars,A,"Af");
   Teko::TpetraHelpers::BlockedTpetraOperator reorderShell(vars,A,"Ar");
 
   Teko::BlockReorderManager brm;
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
   TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ": patern = " << brm.toString());

   TEST_MSG("\n   tBlockedTpetraOperator::test_reorder" << tstr << ":\n");
   TEST_MSG("\n      " << Teuchos::describe(*reorderShell.getThyraOp(), Teuchos::VERB_HIGH)  << std::endl);

   // test the operator against a lot of random vectors
   int numtests = 10;
   ST max = 0.0;
   ST min = 1.0;
   for(int i=0;i<numtests;i++) {
      std::vector<ST> norm(width);
      std::vector<ST> rel(width);
      x.randomize();

      flatShell.apply(x,yf);
      reorderShell.apply(x,yr);

      Tpetra::MultiVector<ST,LO,GO,NT> e(yf,Teuchos::Copy);
      e.update(-1.0,yr,1.0);
      e.norm2(Teuchos::ArrayView<ST>(norm));

      // compute relative error
      yf.norm2(Teuchos::ArrayView<ST>(rel));
      for(int j=0;j<width;j++) {
         max = max>norm[j]/rel[j] ? max : norm[j]/rel[j];
         min = min<norm[j]/rel[j] ? min : norm[j]/rel[j];
      }
   }
   TEST_ASSERT(max>=min,
         "   tBlockedTpetraOperator::test_reorder" << tstr << ": " << toString(status) << ": "
      << "sanity checked - " << max << " >= " << min);
   TEST_ASSERT(max<=tolerance_,
         "   tBlockedTpetraOperator::test_reorder" << tstr << ": "<< toString(status) << ": "
      << "testing tolerance over many matrix vector multiplies ( " << max << " <= "
      << tolerance_ << " )");

   return allPassed;
}

} // end Test namespace
} // end Teko namespace
