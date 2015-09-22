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
#include <Teuchos_TimeMonitor.hpp>

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
#include "Teko_DiagonallyScaledPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_PreconditionerLinearOp.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Test-rig

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

const RCP<Thyra::LinearOpBase<double> > buildSystem(const Epetra_Comm & comm,int size)
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

   return Thyra::nonconstEpetraLinearOp(mat);
}

TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, invfactory_test)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teuchos::ParameterList pl;
   Teuchos::ParameterList & diagList = pl.sublist("DiagScal");
   diagList.set<std::string>("Type","Diagonal Scaling");
   diagList.set<std::string>("Inverse Factory","Amesos");

   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(pl);
   RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("DiagScal"); 
   RCP<Teko::InverseFactory> dirFact = invLib->getInverseFactory("Amesos"); 

   RCP<Thyra::LinearOpBase<double> > A =  buildSystem(Comm,50);

   Teko::LinearOp invA = Teko::buildInverse(*invFact,A);
   Teko::LinearOp invExactA = Teko::buildInverse(*dirFact,A);

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(1e-14);

   const bool result = tester.compare( *invA, *invExactA, Teuchos::ptrFromRef(out));
   if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
   }
   else
      out << "Apply 0: SUCCESS" << std::endl;
}

TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, application_test_row)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   // build linear op tester
   bool result;
   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(1e-14);

   // build operators and factories
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   RCP<Teko::InverseFactory> subsolve = invLib->getInverseFactory("Amesos");
   RCP<Teko::InverseFactory> subsolve_ml = invLib->getInverseFactory("ML");

   RCP<Thyra::LinearOpBase<double> > A =  buildSystem(Comm,50);

   typedef Teko::DiagonallyScaledPreconditionerFactory DSPF;

   RCP<Teko::PreconditionerFactory> dspf = rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve,DSPF::ROW_SCALING));
   RCP<Teko::PreconditionerFactory> dspf_ml = rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve_ml,DSPF::ROW_SCALING));
   RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(dspf,Teuchos::null));
   RCP<Teko::InverseFactory> invFact_ml = rcp(new Teko::PreconditionerInverseFactory(dspf_ml,Teuchos::null));

   // test build inverse capability
   Teko::ModifiableLinearOp invA = Teko::buildInverse(*invFact,A);
   Teko::ModifiableLinearOp invA_ml = Teko::buildInverse(*invFact_ml,A);
   Teko::ModifiableLinearOp invExactA = Teko::buildInverse(*subsolve,A);

   result = tester.compare( *invA, *invExactA, Teuchos::ptrFromRef(out));
   if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
   }
   else
      out << "Apply 0: SUCCESS" << std::endl;

   // test Diagonally scaled rebuild capability
   Teko::rebuildInverse(*invFact,A,invA);
   Teko::rebuildInverse(*invFact_ml,A,invA_ml); // here we tested repeatability using ML

   result = tester.compare( *invA, *invExactA, Teuchos::ptrFromRef(out));
   if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
   }
   else
      out << "Apply 0: SUCCESS" << std::endl;
}

TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, application_test_column)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   // build linear op tester
   bool result;
   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(1e-14);

   // build operators and factories
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   RCP<Teko::InverseFactory> subsolve = invLib->getInverseFactory("Amesos");
   RCP<Teko::InverseFactory> subsolve_ml = invLib->getInverseFactory("ML");

   RCP<Thyra::LinearOpBase<double> > A =  buildSystem(Comm,50);

   RCP<Teko::PreconditionerFactory> dspf = rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve));
   RCP<Teko::PreconditionerFactory> dspf_ml = rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve_ml));
   RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(dspf,Teuchos::null));
   RCP<Teko::InverseFactory> invFact_ml = rcp(new Teko::PreconditionerInverseFactory(dspf_ml,Teuchos::null));

   // test build inverse capability
   Teko::ModifiableLinearOp invA = Teko::buildInverse(*invFact,A);
   Teko::ModifiableLinearOp invA_ml = Teko::buildInverse(*invFact_ml,A);
   Teko::ModifiableLinearOp invExactA = Teko::buildInverse(*subsolve,A);

   result = tester.compare( *invA, *invExactA, Teuchos::ptrFromRef(out));
   if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
   }
   else
      out << "Apply 0: SUCCESS" << std::endl;

   // test Diagonally scaled rebuild capability
   Teko::rebuildInverse(*invFact,A,invA);
   Teko::rebuildInverse(*invFact_ml,A,invA_ml); // here we tested repeatability using ML

   result = tester.compare( *invA, *invExactA, Teuchos::ptrFromRef(out));
   if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
   }
   else
      out << "Apply 0: SUCCESS" << std::endl;
}

TEUCHOS_UNIT_TEST(tDiagonalOperator, replaceValues)
{
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   RCP<Thyra::LinearOpBase<double> > A =  buildSystem(Comm,50);

   Teko::MultiVector diag = Teko::getDiagonal(A,Teko::AbsRowSum);
   Teko::replaceValue(diag,0.0,1.0);
   Teko::LinearOp invDiagOp = Teko::buildInvDiagonal(diag);
}
