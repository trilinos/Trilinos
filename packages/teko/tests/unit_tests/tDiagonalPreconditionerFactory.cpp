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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include <string>
#include <iostream>

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_DiagonalPreconditionerFactory.hpp"
#include "Teko_PreconditionerLinearOp.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

// Test-rig

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

const RCP<Thyra::LinearOpBase<double> > buildSystem(const Epetra_Comm & comm,int size)
{
   Epetra_Map map(-1,size,0,comm);

   RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy,map,0));

   double values[] = { -1.454, 2.0, -1.454}; // modified so lumped matrices don't fail!
   int iTemp[] = {-1,0,1}, indices[3];
   double * vPtr;
   int * iPtr;

   for(int i=0;i<map.NumMyElements();i++) {
      int count = 3;
      int gid = map.GID(i);

      vPtr = values;
      iPtr = indices;

      indices[0] = gid+iTemp[0];
      indices[1] = gid+iTemp[1];
      indices[2] = gid+iTemp[2];
      
      if(gid==0) {
         vPtr = &values[1];
         iPtr = &indices[1];
         count = 2;
      }
      else if(gid==map.MaxAllGID())
         count = 2;

      mat->InsertGlobalValues(gid,count,vPtr,iPtr);
   }

   mat->FillComplete();

   return Thyra::nonconstEpetraLinearOp(mat,"TestOp");
}

RCP<Teuchos::ParameterList> buildLibPL(int blockSize)
{
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   {
      Teuchos::ParameterList & sub = pl->sublist("AbsRowSum");
      sub.set("Type","Explicit Diagonal Preconditioner");
      sub.set("Diagonal Type","AbsRowSum");
   }

   {
      Teuchos::ParameterList & sub = pl->sublist("Diagonal");
      sub.set("Type","Explicit Diagonal Preconditioner");
      sub.set("Diagonal Type","Diagonal");
   }

   {
      Teuchos::ParameterList & sub = pl->sublist("Lumped");
      sub.set("Type","Explicit Diagonal Preconditioner");
      sub.set("Diagonal Type","Lumped");
   }

   {
      Teuchos::ParameterList & sub = pl->sublist("BlkDiag");
      sub.set("Type","Explicit Diagonal Preconditioner");
      sub.set("Diagonal Type","BlkDiag");

      sub.set("contiguous block size",blockSize);
   }

   return pl;
}

TEUCHOS_UNIT_TEST(tDiagonalPreconditionerFactory, diag_inv_test)
{
   using Teuchos::RCP;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   // preconditioners to test
   std::string precName[] = {"AbsRowSum","Diagonal","Lumped"};
   int numPrec = 3;

   RCP<Teuchos::ParameterList> pl = buildLibPL(4);
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*pl); 

   Teko::LinearOp A = buildSystem(Comm,20);
   
   for(int i=0;i<numPrec;i++) {
      RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory(precName[i]);
      Teko::LinearOp idA_fact = Teko::buildInverse(*invFact,A);
      Teko::LinearOp idA_drct = Teko::getInvDiagonalOp(A,Teko::getDiagonalType(precName[i])); 
  
      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
   
      const bool result = tester.compare( *idA_fact, *idA_drct, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE (\"" << precName[i] << "\")" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS (\"" << precName[i] << "\")" << std::endl;

      // test type
      Teko::LinearOp srcOp = Teko::extractOperatorFromPrecOp(idA_fact);
      TEST_ASSERT(Teuchos::rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<double> >(srcOp)!=Teuchos::null);
   }
}

TEUCHOS_UNIT_TEST(tDiagonalPreconditionerFactory, blkdiag_inv_test)
{
   using Teuchos::RCP;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   int numProc = Comm.NumProc();

   // preconditioners to test
   std::string precName = "BlkDiag";

   int myRows = 20;
   int blockSize = 5;
   Teko::LinearOp A = buildSystem(Comm,myRows);

   // sanity check
   TEUCHOS_ASSERT(myRows % blockSize==0);

   {
      RCP<Teuchos::ParameterList> pl = buildLibPL(1); // test diagonal construction
      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*pl); 
   
      RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory(precName);
      Teko::LinearOp idA_fact = Teko::buildInverse(*invFact,A);

      // test type
      Teko::LinearOp srcOp = Teko::extractOperatorFromPrecOp(idA_fact);
      TEST_ASSERT(Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(srcOp)!=Teuchos::null);
   
      Teko::LinearOp idA_drct = Teko::getInvDiagonalOp(A,Teko::Diagonal);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);

      const bool result = tester.compare( *idA_fact, *idA_drct, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE (\"" << precName << "\")" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS (\"" << precName << "\")" << std::endl;
   }

   {
      RCP<Teuchos::ParameterList> pl = buildLibPL(blockSize);
      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*pl); 
   
      RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory(precName);
      Teko::LinearOp idA_fact = Teko::buildInverse(*invFact,A);

      // test type
      Teko::LinearOp srcOp = Teko::extractOperatorFromPrecOp(idA_fact);
      TEST_ASSERT(Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(srcOp)!=Teuchos::null);
      RCP<const Epetra_CrsMatrix> eop = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*srcOp));
      TEST_EQUALITY(eop->NumGlobalNonzeros(),blockSize*blockSize*(numProc*myRows)/blockSize);
      TEST_EQUALITY(eop->NumMyNonzeros(),blockSize*blockSize*myRows/blockSize);
   }
}
