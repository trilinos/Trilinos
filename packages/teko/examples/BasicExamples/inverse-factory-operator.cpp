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
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_InverseFactoryOperator.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // read in parameter list
   Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("simple_example.xml"); 

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
   RCP<Epetra_Vector> b = rcp(ptrb);

   // allocate vectors
   RCP<Epetra_Vector> x = rcp(new Epetra_Vector(A->OperatorDomainMap()));
   x->PutScalar(0.0);

   // Break apart the strided linear system
   /////////////////////////////////////////////////////////

   // Block the linear system using a strided epetra operator
   std::vector<int> vec(2); vec[0] = 2; vec[1] = 1; /*@ \label{lned:define-strided} @*/
   Teuchos::RCP<const Epetra_Operator> strided_A
         = Teuchos::rcp(new Teko::Epetra::StridedEpetraOperator(vec,A));

   // Build the preconditioner /*@ \label{lned:construct-prec} @*/
   /////////////////////////////////////////////////////////

   // build an InverseLibrary and inverse factory
   RCP<Teko::InverseLibrary> invLib
          = Teko::InverseLibrary::buildFromParameterList(*paramList);
   RCP<Teko::InverseFactory> inverse
          = invLib->getInverseFactory("SIMPLE");

   // Create the initial preconditioner, and build it from strided_A
   RCP<Teko::Epetra::InverseFactoryOperator> prec_A
          = rcp(new Teko::Epetra::InverseFactoryOperator(inverse));
   prec_A->initInverse();
   prec_A->rebuildInverseOperator(strided_A.getConst());

   // Build and solve the linear system
   /////////////////////////////////////////////////////////

   // Setup the linear solve: notice A is used directly 
   Epetra_LinearProblem problem(&*A,&*x,&*b); /*@ \label{lned:aztec-solve} @*/

   // build the solver
   AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres);
   solver.SetAztecOption(AZ_precond,AZ_none);
   solver.SetAztecOption(AZ_kspace,1000);
   solver.SetAztecOption(AZ_output,10);
   solver.SetPrecOperator(&*prec_A);

   // solve the linear system
   solver.Iterate(1000,1e-5);

   return 0;
}
