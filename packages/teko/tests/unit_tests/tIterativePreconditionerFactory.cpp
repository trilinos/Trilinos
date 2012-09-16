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
#include "Teko_IterativePreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Test-rig

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

const RCP<const Thyra::LinearOpBase<double> > build2x2(const Epetra_Comm & comm,double a,double b,double c,double d)
{
   RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,comm));

   int indicies[2];
   double row0[2];
   double row1[2];

   indicies[0] = 0;
   indicies[1] = 1;

   // build a CrsMatrix
   RCP<Epetra_CrsMatrix> blk  = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   row0[0] = a; row0[1] = b;  // do a transpose here!
   row1[0] = c; row1[1] = d;
   blk->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   blk->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   blk->FillComplete();

   return Thyra::epetraLinearOp(blk);
}

RCP<Teuchos::ParameterList> buildLibPL(int count,std::string scalingType)
{
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   pl->set("Iterations Count",count);
   pl->set("Preconditioner Type",scalingType);

   return pl;
}

TEUCHOS_UNIT_TEST(tIterativePreconditionerFactory, parameter_list_init)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teko::LinearOp  A = build2x2(Comm,1,2,3,4);

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);

   {
      RCP<Teuchos::ParameterList> pl = buildLibPL(4,"Amesos");
      RCP<Teko::IterativePreconditionerFactory> precFact = rcp(new Teko::IterativePreconditionerFactory());
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact,Teuchos::null));

      try {
         precFact->initializeFromParameterList(*pl);
         out << "Passed correct parameter list" << std::endl;

         Teko::LinearOp prec = Teko::buildInverse(*invFact,A);
      }
      catch(...) {
         success = false; 
         out << "Failed correct parameter list" << std::endl;
      }
   }

   {
      Teuchos::ParameterList pl;
      pl.set("Preconditioner Type","Amesos");

      RCP<Teko::IterativePreconditionerFactory> precFact = rcp(new Teko::IterativePreconditionerFactory());
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact,Teuchos::null));

      try {
         precFact->initializeFromParameterList(pl);
         out << "Passed iteration count" << std::endl;

         Teko::LinearOp prec = Teko::buildInverse(*invFact,A);
      }
      catch(...) {
         out << "Failed iteration count" << std::endl;
      }
   }

   {
      Teuchos::ParameterList pl;
      pl.set("Iteration Count",4);
      pl.set("Precondiioner Type","Amesos");

      RCP<Teko::IterativePreconditionerFactory> precFact = rcp(new Teko::IterativePreconditionerFactory());

      try {
         precFact->initializeFromParameterList(pl);
         success = false;
         out << "Failed preconditioner type" << std::endl;

         // these should not be executed
         RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact,Teuchos::null));
         Teko::LinearOp prec = Teko::buildInverse(*invFact,A);
      }
      catch(const std::exception & exp) {
         out << "Passed preconditioner type" << std::endl;
      }
   }
}

TEUCHOS_UNIT_TEST(tIterativePreconditionerFactory, inverse_test)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teko::LinearOp  A = build2x2(Comm,1,2,3,4);
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");
   Teko::LinearOp iP = Teko::buildInverse(*invFact,A);

   Thyra::LinearOpTester<double> tester;
   tester.dump_all(true);
   tester.show_all_tests(true);

   {
      RCP<Teko::InverseFactory> precOpFact = rcp(new Teko::StaticOpInverseFactory(iP));
      RCP<Teko::IterativePreconditionerFactory> precFact = rcp(new Teko::IterativePreconditionerFactory(9,precOpFact));
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact,Teuchos::null));

      Teko::LinearOp prec = Teko::buildInverse(*invFact,A);

      const bool result = tester.compare( *prec, *iP, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }
}

TEUCHOS_UNIT_TEST(tIterativePreconditionerFactory, constructor_test)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teko::LinearOp  A = build2x2(Comm,1,2,3,4);
   Teko::LinearOp iP = build2x2(Comm,1.0,0.0,0.0,1.0/4.0);
   Teko::LinearOp ImAiP = build2x2(Comm,0.0,-0.5,-3.0,0.0);
   Teko::LinearOp I = Thyra::identity(ImAiP->range());

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);

   {
      RCP<Teko::InverseFactory> precOpFact = rcp(new Teko::StaticOpInverseFactory(iP));
      RCP<Teko::IterativePreconditionerFactory> precFact = rcp(new Teko::IterativePreconditionerFactory(0,precOpFact));
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact,Teuchos::null));

      Teko::LinearOp prec = Teko::buildInverse(*invFact,A);

      const bool result = tester.compare( *prec, *iP, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }

   {
      using Teko::multiply;

      RCP<Teko::InverseFactory> precOpFact = rcp(new Teko::StaticOpInverseFactory(iP));
      RCP<Teko::IterativePreconditionerFactory> precFact = rcp(new Teko::IterativePreconditionerFactory(2,precOpFact));
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact,Teuchos::null));

      Teko::LinearOp prec = Teko::buildInverse(*invFact,A);
      Teko::LinearOp exact = Teko::multiply(iP,Teko::add(I,Teko::multiply(Teko::add(I,ImAiP),ImAiP))); // iP*(I+(I+X)*X)

      const bool result = tester.compare( *prec, *exact, Teuchos::ptrFromRef(out));
      if (!result) {
         out << "Apply 2: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 2: SUCCESS" << std::endl;
   }
}
