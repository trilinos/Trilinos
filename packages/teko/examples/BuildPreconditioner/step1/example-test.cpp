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

// Teuchos includes /*@ \label{lnet:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"

// Epetra includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"

#include "ExamplePreconditionerFactory.cpp"

#include <iostream> /*@ \label{lnet:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

// utility function to construct Epetra operators
RCP<Epetra_CrsMatrix> build2x2(double a,double b,double c,double d,Epetra_Comm & Comm)
{
   RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,Comm));
   RCP<Epetra_CrsMatrix> matrix = rcp(new Epetra_CrsMatrix(Copy,*map,2));

   int indicies[2] = { 0, 1 };
   double values[2];

   // build first row
   values[0] = a;
   values[1] = b;
   matrix->InsertGlobalValues(0,2,values,indicies);

   // build second row
   values[0] = c;
   values[1] = d;
   matrix->InsertGlobalValues(1,2,values,indicies);

   // pack it up
   matrix->FillComplete();

   return matrix;
}

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // build global communicator
   Epetra_MpiComm Comm(MPI_COMM_WORLD);

   // build the sub blocks
   Teko::LinearOp A_00 = Thyra::epetraLinearOp(build2x2(1,2,2,1,Comm));
   Teko::LinearOp A_01 = Thyra::epetraLinearOp(build2x2(0,-1,3,4,Comm));
   Teko::LinearOp A_10 = Thyra::epetraLinearOp(build2x2(1,6,-3,2,Comm));
   Teko::LinearOp A_11 = Thyra::epetraLinearOp(build2x2(2,1,1,2,Comm));

   // build the Epetra operator
   Teuchos::RCP<Teko::Epetra::EpetraOperatorWrapper> A
         = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(Teko::block2x2(A_00,A_01,A_10,A_11)));

   // build the Epetra vector
   Epetra_Vector b(A->OperatorRangeMap());
   Epetra_Vector x(A->OperatorDomainMap());
   x.PutScalar(0.0);

   // build the RHS vector
   b[0] = 1;
   b[1] = 2;
   b[2] = 3;
   b[3] = 4;

   // Build the preconditioner /*@ \label{lnet:construct-prec} @*/
   /////////////////////////////////////////////////////////

   // build an InverseLibrary
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos(); /*@ \label{lnet:define-inv-params} @*/
   
   // build the inverse factory needed by the example preconditioner
   RCP<const Teko::InverseFactory> inverse  /*@ \label{lnet:define-inv-fact} @*/
         = invLib->getInverseFactory("Amesos");

   // build the preconditioner factory
   RCP<Teko::BlockPreconditionerFactory> precFact /*@ \label{lnet:const-prec-fact} @*/
          = rcp(new ExamplePreconditionerFactory(inverse,0.9));

   // using the preconditioner factory construct an Epetra_Operator
   Teko::Epetra::EpetraBlockPreconditioner prec(precFact); /*@ \label{lnet:const-epetra-prec} @*/
   prec.buildPreconditioner(A); // fill epetra preconditioner using the strided operator

   // apply the precondtioner
   prec.ApplyInverse(b,x);

   x.Print(std::cout);

   return 0;
}
