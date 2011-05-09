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

#include "EpetraExt_RowMatrixOut.h"

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactoryOperator.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_StratimikosFactory.hpp"
#include "Teko_StridedEpetraOperator.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

#include "Teuchos_AbstractFactoryStd.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

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

   return mat;
}

const RCP<Epetra_Operator> buildStridedSystem(const Epetra_Comm & comm,int size)
{
   Epetra_Map map(2*size,0,comm);

   RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy,map,0));

   int numUnks = 2;
   double valuesA[] = { -1.0, 2.0, 7.0, -1.0 };
   double valuesB[] = { -1.0, 2.0, 9.0, -1.0 };
   int iTempA[] = {-numUnks,0, 1,numUnks}, indices[4];
   int iTempB[] = {-numUnks,0,-1,numUnks};
   double * vPtr;
   int * iPtr;

   for(int i=0;i<map.NumMyElements()/numUnks;i++) {
      int count = 4;
      int gidA = map.GID(2*i);
      int gidB = gidA+1;

      for(int n=0;n<numUnks;n++) {
         int * iTemp = (n==0) ? iTempA : iTempB;
         int gid = (n==0) ? gidA : gidB;

         indices[0] = gid+iTemp[0];
         indices[1] = gid+iTemp[1];
         indices[2] = gid+iTemp[2];
         indices[3] = gid+iTemp[3];

         vPtr = (n==0) ? valuesA : valuesB;
         iPtr = indices;
      
         if(gid<numUnks) {
            vPtr++;
            iPtr = &indices[1];
            count = 3;
         }
         else if(gid>map.MaxAllGID()-numUnks)
            count = 3;

         mat->InsertGlobalValues(gid,count,vPtr,iPtr);
      }
   }

   mat->FillComplete();

   EpetraExt::RowMatrixToMatrixMarketFile("strided.mm",*mat);

   return mat;
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_Defaults) 
{
   using Teuchos::RCP;
   using Teuchos::ParameterList;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm comm;
   #endif

   // build epetra operator
   RCP<Epetra_Operator> eA = buildSystem(comm,5);
   RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

   // build stratimikos factory, adding Teko's version
   Stratimikos::DefaultLinearSolverBuilder stratFactory;
   stratFactory.setPreconditioningStrategyFactory(
         Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,Teko::StratimikosFactory>(),
         "Teko");
   RCP<const ParameterList> validParams = stratFactory.getValidParameters();
   stratFactory.setParameterList(Teuchos::rcp(new Teuchos::ParameterList(*validParams)));

   // print out Teko's parameter list and fail if it doesn't exist!
   TEST_NOTHROW(validParams->sublist("Preconditioner Types").sublist("Teko").print(out,
         ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true)));

   // build teko preconditioner factory
   RCP<Thyra::PreconditionerFactoryBase<double> > precFactory
         = stratFactory.createPreconditioningStrategy("Teko");

   // make sure factory is built
   TEST_ASSERT(precFactory!=Teuchos::null);

   // build preconditioner
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,tA);
   TEST_ASSERT(prec!=Teuchos::null);

   // build an operator to test against
   RCP<const Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
   RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");
   Teko::LinearOp testOp = Teko::buildInverse(*invFact,tA);

   Teko::LinearOp precOp = prec->getUnspecifiedPrecOp();
   TEST_ASSERT(precOp!=Teuchos::null);

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(0);
   TEST_ASSERT(tester.compare(*precOp,*testOp,Teuchos::ptrFromRef(out)));
}


TEUCHOS_UNIT_TEST(tStratimikosFactory, test_BlockGaussSeidel) 
{
   using Teuchos::RCP;
   using Teuchos::ParameterList;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm comm;
   #endif

   // build epetra operator
   RCP<Epetra_Operator> eA = buildStridedSystem(comm,5);
   RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

   // build stratimikos factory, adding Teko's version
   Stratimikos::DefaultLinearSolverBuilder stratFactory;
   stratFactory.setPreconditioningStrategyFactory(
         Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,Teko::StratimikosFactory>(),
         "Teko");
   RCP<ParameterList> params = Teuchos::rcp(new ParameterList(*stratFactory.getValidParameters()));
   ParameterList & tekoList = params->sublist("Preconditioner Types").sublist("Teko");
   tekoList.set("Write Block Operator", false);
   tekoList.set("Test Block Operator", false);
   tekoList.set("Strided Blocking","1 1");
   tekoList.set("Inverse Type","BGS");
   ParameterList & ifl = tekoList.sublist("Inverse Factory Library");
   ifl.sublist("BGS").set("Type","Block Gauss-Seidel");
   ifl.sublist("BGS").set("Inverse Type","Amesos");

   // RCP<Thyra::PreconditionerFactoryBase<double> > precFactory
   //       = stratFactory.createPreconditioningStrategy("Teko");

   // build operator to test against
   Teko::LinearOp testOp;
   {
      Teuchos::ParameterList iflCopy(ifl);
      RCP<Epetra_Operator> strided_eA = Teuchos::rcp(new Teko::Epetra::StridedEpetraOperator(2,eA));
      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(iflCopy);
      RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("BGS");
      RCP<Teko::Epetra::InverseFactoryOperator> invFactOp = Teuchos::rcp(new Teko::Epetra::InverseFactoryOperator(invFact));
      invFactOp->buildInverseOperator(strided_eA);

      testOp = Thyra::epetraLinearOp(invFactOp,Thyra::NOTRANS,Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);
   }

   stratFactory.setParameterList(params);
   RCP<Thyra::PreconditionerFactoryBase<double> > precFactory
         = stratFactory.createPreconditioningStrategy("Teko");

   // build teko preconditioner factory
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,tA);
   Teko::LinearOp precOp = prec->getUnspecifiedPrecOp();
   TEST_ASSERT(precOp!=Teuchos::null);

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   tester.set_all_error_tol(0);
   TEST_ASSERT(tester.compare(*precOp,*testOp,Teuchos::ptrFromRef(out)));
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_RelatedFunctions) 
{
   using Teuchos::RCP;
   using Teuchos::ParameterList;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm comm;
   #endif

   // build epetra operator
   RCP<Epetra_Operator> eA = buildStridedSystem(comm,5);
   RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

   {
      // build stratimikos factory, adding Teko's version
      Stratimikos::DefaultLinearSolverBuilder stratFactory;
      
      Teko::addTekoToStratimikosBuilder(stratFactory);
      TEST_THROW(Teko::addTekoToStratimikosBuilder(stratFactory),std::logic_error);
      Teko::addTekoToStratimikosBuilder(stratFactory,"Teko-2");

      TEST_NOTHROW(stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko"));
      TEST_NOTHROW(stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko-2"));
   }

   {
      Teuchos::RCP<Teko::RequestHandler> rh = Teuchos::rcp(new Teko::RequestHandler);

      // build stratimikos factory, adding Teko's version
      Stratimikos::DefaultLinearSolverBuilder stratFactory;
      
      Teko::addTekoToStratimikosBuilder(stratFactory,rh);
      TEST_THROW(Teko::addTekoToStratimikosBuilder(stratFactory,rh),std::logic_error);
      Teko::addTekoToStratimikosBuilder(stratFactory,rh,"Teko-2");

      TEST_NOTHROW(stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko"));
      TEST_NOTHROW(stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko-2"));

      RCP<ParameterList> params = Teuchos::rcp(new ParameterList(*stratFactory.getValidParameters()));
      ParameterList & tekoList = params->sublist("Preconditioner Types").sublist("Teko");
      tekoList.set("Write Block Operator", false);
      tekoList.set("Test Block Operator", false);
      tekoList.set("Strided Blocking","1 1");
      tekoList.set("Inverse Type","BGS");
      ParameterList & ifl = tekoList.sublist("Inverse Factory Library");
      ifl.sublist("BGS").set("Type","Block Gauss-Seidel");
      ifl.sublist("BGS").set("Inverse Type","Amesos");
      stratFactory.setParameterList(params);

      RCP<Teko::StratimikosFactory> precFactory
            = Teuchos::rcp_dynamic_cast<Teko::StratimikosFactory>(stratFactory.createPreconditioningStrategy("Teko-2"));
      TEST_EQUALITY(precFactory->getRequestHandler(),rh);
   }
}
