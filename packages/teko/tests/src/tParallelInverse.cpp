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

#include "tParallelInverse.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedEpetraOperator.hpp"

#include "Thyra_EpetraLinearOp.hpp"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::epetraLinearOp;

void tParallelInverse::initializeTest()
{
   tolerance_ = 1.0e-7;
}

void tParallelInverse::loadMatrix()
{
   // Read in the matrix, store pointer as an RCP
   Epetra_CrsMatrix * ptrA = 0;
   TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToCrsMatrix("data/nsjac_1.mm",*GetComm(),ptrA));
   F_ = Thyra::epetraLinearOp(rcp(ptrA));
}

void tParallelInverse::loadStridedMatrix()
{
   // Read in the matrix, store pointer as an RCP
   Epetra_CrsMatrix * ptrA = 0;
   TEST_FOR_EXCEPT(EpetraExt::MatrixMarketFileToCrsMatrix("data/nsjac.mm",*GetComm(),ptrA));
   RCP<Epetra_CrsMatrix> A = rcp(ptrA);

   // Block the linear system using a strided epetra operator
   std::vector<int> vec(2); vec[0] = 1; vec[1] = 2; /*@ \label{lned:define-strided} @*/
   Teko::Epetra::StridedEpetraOperator sA(vec,A);

   // get 0,0 block
   F_ = Thyra::epetraLinearOp(sA.GetBlock(0,0));
}

int tParallelInverse::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tParallelInverse";


   status = test_inverse(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"withmassStable\" ... PASSED","   \"withmassStable\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_stridedInverse(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"nomassStable\" ... PASSED","   \"nomassStable\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tParallelInverse...PASSED","tParallelInverse...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tParallelInverse...FAILED");
   }

   return failcount;
}

bool tParallelInverse::test_inverse(int verbosity,std::ostream & os)
{
   // bool status = false;
   bool allPassed = true;

   loadMatrix();

   // build an InverseLibrary
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

   // build the inverse factory needed by the example preconditioner
   RCP<const Teko::InverseFactory> invFact  = invLib->getInverseFactory("Amesos");

   Teko::LinearOp inv = invFact->buildInverse(F_);

   return allPassed;
}

bool tParallelInverse::test_stridedInverse(int verbosity,std::ostream & os)
{
   // bool status = false;
   bool allPassed = true;

   loadStridedMatrix();

   // build an InverseLibrary
   TEST_MSG("\n   Building inverse library");
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

   // build the inverse factory needed by the example preconditioner
   TEST_MSG("\n   Building inverse factory");
   RCP<const Teko::InverseFactory> invFact  = invLib->getInverseFactory("Amesos");

   TEST_MSG("\n   Building inverse");
   Teko::LinearOp inv = invFact->buildInverse(F_);

   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
