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

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactoryOperator.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#define SS_ECHO(ops) { std::stringstream ss; ss << ops; ECHO(ss.str()) };

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

///////////////////////////////////////////////////////////

const RCP<Epetra_Operator> buildSystem(const Epetra_Comm & comm,int size)
{
   Epetra_Map map(size,0,comm);

   RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy,map,0));

   double values[] = { -1.0, 2.0, -1.0};
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

   // return Thyra::nonconstEpetraLinearOp(mat);
   return mat;
}

TEUCHOS_UNIT_TEST(tInverseFactoryOperator, test_Direct_Solve) 
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm comm;
   #endif

   Teuchos::RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   Teuchos::RCP<Teko::InverseFactory> invFactory
         = invLib->getInverseFactory("Amesos");
  
   Teuchos::RCP<Epetra_Operator> eA = buildSystem(comm,50);
   Teko::LinearOp A = Thyra::epetraLinearOp(eA);
   Teko::ModifiableLinearOp invA = Teko::buildInverse(*invFactory,A);

   Teko::Epetra::InverseFactoryOperator invFactOp(invFactory);
   invFactOp.buildInverseOperator(eA);

   {
      // because InverseFactoryOperator is a "Preconditioner" then need to
      // call Epetra_Operator::ApplyInverse
      Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp),
                                   Thyra::NOTRANS,Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
      tester.set_all_error_tol(1e-14);
   
      const bool result = tester.compare( *invA, *testInvA, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }

   invFactOp.rebuildInverseOperator(eA);
   {
      // because InverseFactoryOperator is a "Preconditioner" then need to
      // call Epetra_Operator::ApplyInverse
      Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp),
                                   Thyra::NOTRANS,Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
      tester.set_all_error_tol(1e-14);
   
      const bool result = tester.compare( *invA, *testInvA, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }
}

TEUCHOS_UNIT_TEST(tInverseFactoryOperator, test_Block_Solve) 
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm comm;
   #endif

   Teuchos::RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   Teuchos::RCP<Teko::InverseFactory> amesosFactory
         = invLib->getInverseFactory("Amesos");

   Teuchos::RCP<Epetra_Operator> eA00 = buildSystem(comm,50);
   Teko::LinearOp A_00 = Thyra::epetraLinearOp(eA00);
   Teko::ModifiableLinearOp invA_00 = Teko::buildInverse(*amesosFactory,A_00);

   Teko::LinearOp A = Thyra::block2x2<double>(A_00,Teuchos::null,Teuchos::null,A_00);
   Teko::LinearOp invA = Thyra::block2x2<double>(invA_00,Teuchos::null,Teuchos::null,invA_00); 
   Teuchos::RCP<Epetra_Operator> eInvA = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(invA));
   Teko::LinearOp cmpInvA = Thyra::epetraLinearOp(eInvA);

   Teuchos::RCP<Teko::PreconditionerFactory> jacFact = Teuchos::rcp(new Teko::JacobiPreconditionerFactory(invA_00,invA_00));
   Teuchos::RCP<Teko::InverseFactory> invFactory = Teuchos::rcp(new Teko::PreconditionerInverseFactory(jacFact,Teuchos::null));
   Teuchos::RCP<Epetra_Operator> eA = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(A));

   Teko::Epetra::InverseFactoryOperator invFactOp(invFactory);
   invFactOp.buildInverseOperator(eA);

   {
      // because InverseFactoryOperator is a "Preconditioner" then need to
      // call Epetra_Operator::ApplyInverse
      Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp),
                                   Thyra::NOTRANS,Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
      tester.set_all_error_tol(1e-14);
   
      const bool result = tester.compare( *cmpInvA, *testInvA, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }

   invFactOp.rebuildInverseOperator(eA);
   {
      // because InverseFactoryOperator is a "Preconditioner" then need to
      // call Epetra_Operator::ApplyInverse
      Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp),
                                   Thyra::NOTRANS,Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
      tester.set_all_error_tol(1e-14);
   
      const bool result = tester.compare( *cmpInvA, *testInvA, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }
}
