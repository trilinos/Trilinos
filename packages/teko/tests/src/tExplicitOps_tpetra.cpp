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

#include "tExplicitOps_tpetra.hpp"

#include <string>

// Tpetra includes
#include "Tpetra_Export.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

void tExplicitOps_tpetra::initializeTest()
{
   const Epetra_Comm & comm_epetra = *GetComm();
   const RCP<const Teuchos::Comm<int> > comm_tpetra = GetComm_tpetra();

   tolerance_ = 1.0e-4;

   int nx = 39; // essentially random values
   int ny = 53;

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm_epetra,false);
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   Epetra_CrsMatrix & epetraF = FGallery.GetMatrixRef();
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > tpetraF = Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(rcpFromRef(epetraF),comm_tpetra);
   F_ = Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraF->getDomainMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraF->getRangeMap()),tpetraF);

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery GGallery("laplace_2d",comm_epetra,false);
   GGallery.Set("nx",nx);
   GGallery.Set("ny",ny);
   Epetra_CrsMatrix & epetraG = GGallery.GetMatrixRef();
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > tpetraG = Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(rcpFromRef(epetraG),comm_tpetra);
   G_ = Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraG->getDomainMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tpetraG->getRangeMap()),tpetraG);

   RCP<Tpetra::Vector<ST,LO,GO,NT> > v = rcp(new Tpetra::Vector<ST,LO,GO,NT> (tpetraF->getRangeMap()));
   v->randomize();
   RCP<Thyra::VectorBase<ST> > tV = Thyra::createVector<ST,LO,GO,NT>(v,Thyra::createVectorSpace<ST,LO,GO,NT>(tpetraF->getRowMap())); 
   D_ = Thyra::diagonal(tV);
}

int tExplicitOps_tpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tExplicitOps_tpetra";

   status = test_mult_diagScaleMatProd(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"mult_diagScaleMatProd\" ... PASSED","   \"mult_diagScaleMatProd\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_mult_diagScaling(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"mult_diagScaling\" ... PASSED","   \"mult_diagScaling\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_add(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"add\" ... PASSED","   \"add\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_mult_modScaleMatProd(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"mult_modScaleMatProd\" ... PASSED","   \"mult_modScaleMatProd\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_add_mod(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"add_mod\" ... PASSED","   \"add\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tExplicitOps_tpetra...PASSED","tExplicitOps_tpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tExplicitOps_tpetra...FAILED");
   }

   return failcount;
}

bool tExplicitOps_tpetra::test_mult_diagScaleMatProd(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<ST> > thyOp;
   Teko::LinearOp expOp;

   thyOp = Thyra::multiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_diagScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   thyOp = Teko::multiply(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),Teko::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_diagScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps_tpetra::test_mult_diagScaling(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<ST> > thyOp;
   Teko::LinearOp expOp;

   thyOp = Teko::multiply(Teko::scale(-4.0,F_),D_);
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_diagScaleMatProd "
             << ": Testing diagonal scaling");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   thyOp = Teko::multiply(D_,Teko::scale(-9.0,F_));
   expOp = Teko::explicitMultiply(D_,Teko::scale(-9.0,F_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_diagScaleMatProd "
             << ": Testing diagonal scaling");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps_tpetra::test_mult_modScaleMatProd(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   Teko::LinearOp thyOp;
   Teko::ModifiableLinearOp expOp;

   thyOp = Teko::multiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_),expOp);

   RCP<const Thyra::TpetraLinearOp<ST,LO,GO,NT> > tOp1 = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST,LO,GO,NT> >(expOp,true);
   RCP<const Tpetra::Operator<ST,LO,GO,NT> > eop1 = tOp1->getConstTpetraOperator();

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_modScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   RCP<Thyra::TpetraLinearOp<ST,LO,GO,NT> > tF = Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST,LO,GO,NT> >(F_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > crsF = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST,LO,GO,NT> >(tF->getTpetraOperator(),true);
   crsF->resumeFill();
   crsF->scale(5.0);
   crsF->fillComplete();
   RCP<Thyra::TpetraLinearOp<ST,LO,GO,NT> > tG = Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST,LO,GO,NT> >(G_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > crsG = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST,LO,GO,NT> >(tG->getTpetraOperator(),true);
   crsG->resumeFill();
   crsG->scale(2.0);
   crsG->fillComplete();

   // do some random violence (oh my brothers) to one row
   size_t numEntries = crsF->getNumEntriesInLocalRow (3);
   auto indices1 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_local_inds_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   auto values1 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_values_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   crsF->getLocalRowCopy(3,indices1,values1,numEntries);
   for(size_t i=0;i<numEntries;i++) values1(i) *= values1(i)*ST(i+1)*0.92;
   crsF->replaceLocalValues(3,indices1,values1);

   // do some random violence (oh my brothers) to one row
   numEntries = crsF->getNumEntriesInLocalRow (7);
   auto indices2 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_local_inds_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   auto values2 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_values_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   crsF->getLocalRowCopy(7,indices2,values2,numEntries);
   for(size_t i=0;i<numEntries;i++) values2(i) *= values2(i)*ST(i+1)*0.92;
   crsF->replaceLocalValues(7,indices2,values2);

   // perform the next test
   thyOp = Teko::multiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_),expOp);

   RCP<const Thyra::TpetraLinearOp<ST,LO,GO,NT> > tOp2 = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST,LO,GO,NT> >(expOp,true);
   RCP<const Tpetra::Operator<ST,LO,GO,NT> > eop2 = tOp2->getConstTpetraOperator();

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_modScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps_tpetra::test_add(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<ST> > thyOp;
   Teko::LinearOp expOp;

   thyOp = Teko::add(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitAdd(Teko::scale(-4.0,F_),Teko::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_add "
             << ": Testing explicit add");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps_tpetra::test_add_mod(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<ST> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<ST> > thyOp;
   Teko::ModifiableLinearOp expOp;

   thyOp = Teko::add(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitAdd(Teko::scale(-4.0,F_),Teko::adjoint(G_),expOp);

   RCP<const Thyra::TpetraLinearOp<ST,LO,GO,NT> > tOp1 = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST,LO,GO,NT> >(expOp,true);
   RCP<const Tpetra::Operator<ST,LO,GO,NT> > eop1 = tOp1->getConstTpetraOperator();

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_add_mod"
             << ": Testing explicit add");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   RCP<Thyra::TpetraLinearOp<ST,LO,GO,NT> > tF = Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST,LO,GO,NT> >(F_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > crsF = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST,LO,GO,NT> >(tF->getTpetraOperator(),true);
   crsF->resumeFill();
   crsF->scale(5.0);
   crsF->fillComplete();
   RCP<Thyra::TpetraLinearOp<ST,LO,GO,NT> > tG = Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST,LO,GO,NT> >(G_);
   RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > crsG = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST,LO,GO,NT> >(tG->getTpetraOperator(),true);
   crsG->resumeFill();
   crsG->scale(2.0);
   crsG->fillComplete();

   // do some random violence (oh my brothers) to one row
   size_t numEntries = crsF->getNumEntriesInLocalRow (3);
   auto indices1 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_local_inds_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   auto values1 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_values_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   crsF->getLocalRowCopy(3,indices1,values1,numEntries);
   for(size_t i=0;i<numEntries;i++) values1(i) *= values1(i)*ST(i+1)*0.92;
   crsF->replaceLocalValues(3,indices1,values1);

   // do some random violence (oh my brothers) to one row
   numEntries = crsF->getNumEntriesInLocalRow (7);
   auto indices2 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_local_inds_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   auto values2 = typename Tpetra::CrsMatrix<ST,LO,GO,NT>::nonconst_values_host_view_type(Kokkos::ViewAllocateWithoutInitializing("rowIndices"),numEntries);
   crsF->getLocalRowCopy(7,indices2,values2,numEntries);
   for(size_t i=0;i<numEntries;i++) values2(i) *= values2(i)*ST(i+1)*0.92;
   crsF->replaceLocalValues(7,indices2,values2);

   // perform the next test
   thyOp = Teko::add(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitAdd(Teko::scale(-4.0,F_),Teko::adjoint(G_),expOp);

   RCP<const Thyra::TpetraLinearOp<ST,LO,GO,NT> > tOp2 = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST,LO,GO,NT> >(expOp,true);
   RCP<const Tpetra::Operator<ST,LO,GO,NT> > eop2 = tOp2->getConstTpetraOperator();

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps_tpetra::test_add_mod"
             << ": Testing matrix addition");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
