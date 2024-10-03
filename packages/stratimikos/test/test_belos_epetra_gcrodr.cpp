// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <string>
#include <iostream>

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
// #include "EpetraExt_RowMatrixOut.h"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveHelpers.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

Teuchos::RCP<Epetra_CrsMatrix> buildMatrix(int nx, Epetra_Comm & comm)
{
   Epetra_Map map(nx*comm.NumProc(),0,comm);
   Teuchos::RCP<Epetra_CrsMatrix> mat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,map,3));

   int offsets[3] = {-1, 0, 1 };
   double values[3] = { -1, 2, -1};
   int maxGid = map.MaxAllGID();
   for(int lid=0;lid<nx;lid++) {
      int gid = mat->GRID(lid);
      int numEntries = 3, offset = 0;
      int indices[3] = { gid+offsets[0],
                         gid+offsets[1],
                         gid+offsets[2] };
      if(gid==0) { // left end point
         numEntries = 2;
         offset = 1;
      }            // right end point
      else if(gid==maxGid)
         numEntries = 2;

      // insert rows
      mat->InsertGlobalValues(gid,numEntries,values+offset,indices+offset);
   }

   mat->FillComplete();
   return mat;
}

TEUCHOS_UNIT_TEST(belos_gcrodr, multiple_solves)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   // build and allocate linear system
   Teuchos::RCP<Epetra_CrsMatrix> mat = buildMatrix(100,Comm);
   Teuchos::RCP<Epetra_Vector> x0 = rcp(new Epetra_Vector(mat->OperatorDomainMap()));
   Teuchos::RCP<Epetra_Vector> x1 = rcp(new Epetra_Vector(mat->OperatorDomainMap()));
   Teuchos::RCP<Epetra_Vector> b = rcp(new Epetra_Vector(mat->OperatorRangeMap()));

   x0->Random();
   x1->Random();
   b->PutScalar(0.0);

   // sanity check
   // EpetraExt::RowMatrixToMatrixMarketFile("mat_output.mm",*mat);

   // build Thyra wrappers
   RCP<const Thyra::LinearOpBase<double> >
      tA = Thyra::epetraLinearOp( mat );
   RCP<Thyra::VectorBase<double> >
      tx0 = Thyra::create_Vector( x0, tA->domain() );
   RCP<Thyra::VectorBase<double> >
      tx1 = Thyra::create_Vector( x1, tA->domain() );
   RCP<const Thyra::VectorBase<double> >
      tb = Thyra::create_Vector( b, tA->range() );

   // now comes Stratimikos
   RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("BelosGCRODRTest.xml");
   Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
   linearSolverBuilder.setParameterList(paramList);
 
   // Create a linear solver factory given information read from the
   // parameter list.
   RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
         linearSolverBuilder.createLinearSolveStrategy("");

   // Create a linear solver based on the forward operator A
   RCP<Thyra::LinearOpWithSolveBase<double> > lows =
         Thyra::linearOpWithSolve(*lowsFactory, tA);

   // Solve the linear system 
   Thyra::SolveStatus<double> status; 
   status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *tb, tx0.ptr());
   status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *tb, tx1.ptr());
}

TEUCHOS_UNIT_TEST(belos_gcrodr, 2x2_multiple_solves)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   // build and allocate linear system
   Teuchos::RCP<Epetra_CrsMatrix> mat = buildMatrix(100,Comm);
   Teuchos::RCP<Epetra_Vector> b = rcp(new Epetra_Vector(mat->OperatorRangeMap()));

   b->PutScalar(0.0);

   // sanity check
   // EpetraExt::RowMatrixToMatrixMarketFile("mat_output.mm",*mat);

   // build Thyra wrappers
   RCP<const Thyra::LinearOpBase<double> > tA;
   RCP<const Thyra::VectorBase<double> > tb;
   {
      // build blocked linear Op
      RCP<const Thyra::LinearOpBase<double> > tA_sub 
            = Thyra::epetraLinearOp( mat );
      RCP<const Thyra::LinearOpBase<double> > zero 
            = Thyra::zero(tA_sub->range(),tA_sub->domain());
      
      tA = Thyra::block2x2(tA_sub,zero,zero,tA_sub);

      // build blocked vector
      RCP<const Thyra::VectorBase<double> > tb_sub 
            = Thyra::create_Vector( b, tA_sub->range() );

      RCP<Thyra::VectorBase<double> > tb_m = Thyra::createMember(tA->range());
      Thyra::randomize(-1.0,1.0,tb_m.ptr());

      tb = tb_m;
   }
   RCP<Thyra::VectorBase<double> > tx0;
   RCP<Thyra::VectorBase<double> > tx1;
   {
      tx0 = Thyra::createMember(tA->domain());
      tx1 = Thyra::createMember(tA->domain());
     
      Thyra::randomize(-1.0,1.0,tx0.ptr());
      Thyra::randomize(-1.0,1.0,tx1.ptr());
   }

   // now comes Stratimikos
   RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("BelosGCRODRTest.xml");
   Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
   linearSolverBuilder.setParameterList(paramList);
 
   // Create a linear solver factory given information read from the
   // parameter list.
   RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
         linearSolverBuilder.createLinearSolveStrategy("");

   // Create a linear solver based on the forward operator A
   RCP<Thyra::LinearOpWithSolveBase<double> > lows =
         Thyra::linearOpWithSolve(*lowsFactory, tA);

   // Solve the linear system 
   Thyra::SolveStatus<double> status; 
   status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *tb, tx0.ptr());
   status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *tb, tx1.ptr());
}
