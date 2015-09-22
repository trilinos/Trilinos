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

#include "tLumping.hpp"

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
   Teko_TEST_MSG(stdstrm,1,"   \"lumping\" ... PASSED","   \"lumping\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_invLumping(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"invLumping\" ... PASSED","   \"invLumping\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tLumping...PASSED","tLumping...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tLumping...FAILED");
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

   Teko::LinearOp pA = Thyra::epetraLinearOp(rcpFromRef(A));
   Teko::LinearOp pB = Thyra::epetraLinearOp(rcpFromRef(B));
   Teko::LinearOp lumpedA = getLumpedMatrix(pA);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *pB, *lumpedA, Teuchos::ptrFromRef(fos) );
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

   Teko::LinearOp pA = Thyra::epetraLinearOp(rcpFromRef(A));
   Teko::LinearOp pB = Thyra::epetraLinearOp(rcpFromRef(B));
   Teko::LinearOp lumpedA = getInvLumpedMatrix(pA);

   {
      std::stringstream ss;
      Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");
      const bool result = tester.compare( *pB, *lumpedA, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tLummping::test_invLumping "
             << ": Testing basic inv lumping functionality");
      if(not result || verbosity>=10)
         os << ss.str();
   }

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
