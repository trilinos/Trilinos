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

// Teuchos includes /*@ \label{lned:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Epetra includes
#ifdef HAVE_MPI
   #include "mpi.h"
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_StratimikosFactory.hpp"

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_VectorBase.hpp"

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // read in parameter list
   Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("strat_example.xml"); 

   // build global communicator
#ifdef HAVE_MPI
   Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
   Epetra_SerialComm Comm;
#endif

   // Read in the matrix, store pointer as an RCP
   Epetra_CrsMatrix * ptrA = 0;
   EpetraExt::MatrixMarketFileToCrsMatrix("../data/nsjac_test.mm",Comm,ptrA);
   RCP<Epetra_CrsMatrix> A = rcp(ptrA);

   // read in the RHS vector
   Epetra_Vector * ptrb = 0;
   EpetraExt::MatrixMarketFileToVector("../data/nsrhs_test.mm",A->OperatorRangeMap(),ptrb);
   RCP<const Epetra_Vector> b = rcp(ptrb);

   // allocate vectors
   RCP<Epetra_Vector> x = rcp(new Epetra_Vector(A->OperatorDomainMap()));
   x->PutScalar(0.0);

   // Build Thyra linear algebra objects
   RCP<const Thyra::LinearOpBase<double> > th_A = Thyra::epetraLinearOp(A);
   RCP<const Thyra::VectorBase<double> > th_b = Thyra::create_Vector(b,th_A->range());
   RCP<Thyra::VectorBase<double> > th_x = Thyra::create_Vector(x,th_A->domain());

   // Build stratimikos solver
   /////////////////////////////////////////////////////////

   Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

   Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
   linearSolverBuilder.setParameterList(paramList);

   RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory 
       = Thyra::createLinearSolveStrategy(linearSolverBuilder);

   Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > th_invA 
      = Thyra::linearOpWithSolve(*lowsFactory, th_A);

   Thyra::assign(th_x.ptr(), 0.0);
   Thyra::SolveStatus<double> status 
       = Thyra::solve<double>(*th_invA, Thyra::NOTRANS, *th_b, th_x.ptr());

   return 0;
}
