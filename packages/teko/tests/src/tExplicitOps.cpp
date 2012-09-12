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

#include "tExplicitOps.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"

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

namespace Teko {
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
      Teko_TEST_MSG(failstrm,0,"tExplicitOps...PASSED","tExplicitOps...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tExplicitOps...FAILED");
   }

   return failcount;
}

bool tExplicitOps::test_mult_diagScaleMatProd(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<double> > thyOp;
   Teko::LinearOp expOp;

   thyOp = Thyra::multiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_diagScaleMatProd "
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
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<double> > thyOp;
   Teko::LinearOp expOp;

   thyOp = Teko::multiply(Teko::scale(-4.0,F_),D_);
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_diagScaleMatProd "
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
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   Teko::LinearOp thyOp;
   Teko::ModifiableLinearOp expOp;

   thyOp = Teko::multiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_),expOp);

   RCP<const Epetra_Operator> eop1 = Thyra::get_Epetra_Operator(*expOp);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
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
   thyOp = Teko::multiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_));
   expOp = Teko::explicitMultiply(Teko::scale(-4.0,F_),D_,Teko::adjoint(G_),expOp);

   RCP<const Epetra_Operator> eop2 = Thyra::get_Epetra_Operator(*expOp);
   TEUCHOS_ASSERT(eop1==eop2);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_modScaleMatProd "
             << ": Testing triple matrix product");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps::test_add(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<double> > thyOp;
   Teko::LinearOp expOp;

   thyOp = Teko::add(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitAdd(Teko::scale(-4.0,F_),Teko::adjoint(G_));

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_add "
             << ": Testing explicit add");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

bool tExplicitOps::test_add_mod(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   Thyra::LinearOpTester<double> tester;
   tester.set_all_error_tol(1e-10);
   tester.show_all_tests(true);

   RCP<const Thyra::LinearOpBase<double> > thyOp;
   Teko::ModifiableLinearOp expOp;

   thyOp = Teko::add(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitAdd(Teko::scale(-4.0,F_),Teko::adjoint(G_),expOp);

   RCP<const Epetra_Operator> eop1 = Thyra::get_Epetra_Operator(*expOp);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_add_mod"
             << ": Testing explicit add");
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
   thyOp = Teko::add(Teko::scale(-4.0,F_),Teko::adjoint(G_));
   expOp = Teko::explicitAdd(Teko::scale(-4.0,F_),Teko::adjoint(G_),expOp);

   RCP<const Epetra_Operator> eop2 = Thyra::get_Epetra_Operator(*expOp);
   TEUCHOS_ASSERT(eop1==eop2);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *thyOp, *expOp, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tExplicitOps::test_add_mod"
             << ": Testing matrix addition");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
