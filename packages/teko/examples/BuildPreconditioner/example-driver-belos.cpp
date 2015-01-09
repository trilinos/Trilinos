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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

// Epetra includes
#include "mpi.h"

#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "MatrixMarket_Tpetra.hpp"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_Ifpack2PreconditionerFactory.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosThyraAdapter.hpp" // Requires Stratimikos...

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

RCP<Teuchos::ParameterList> buildLibPL();

int main(int argc,char * argv[])
{
   typedef double                            ST;
   typedef Thyra::MultiVectorBase<ST>        MV;
   typedef Thyra::LinearOpBase<ST>           OP;

   typedef KokkosClassic::DefaultNode::DefaultNodeType NT;
   typedef Tpetra::CrsMatrix<ST,int,int> TP_Crs;
   typedef Tpetra::Operator<ST,int,int> TP_Op;
   typedef Tpetra::Vector<ST,int,int> TP_Vec;

   typedef Thyra::PreconditionerFactoryBase<double> Base;
   typedef Thyra::Ifpack2PreconditionerFactory<TP_Crs> Impl;

   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // tell Stratimikos => Teko about Ifpack2
   RCP<Stratimikos::DefaultLinearSolverBuilder> linearSolverBuilder = Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder);
   linearSolverBuilder->setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");

   // read in the CRS matrix
   RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
   RCP<TP_Op> Mat = Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile("../data/nsjac_test.mm",
                                                                         Teuchos::DefaultComm<int>::getComm(),
                                                                         node);

   // Allocate some right handside vectors
   RCP<TP_Vec> x0_tp = rcp(new TP_Vec(Mat->getDomainMap()));
   RCP<TP_Vec> x1_tp = rcp(new TP_Vec(Mat->getDomainMap()));
   RCP<TP_Vec> b0_tp = rcp(new TP_Vec(Mat->getRangeMap()));
   RCP<TP_Vec> b1_tp = rcp(new TP_Vec(Mat->getRangeMap()));
   b0_tp->randomize();
   b1_tp->randomize();

   RCP<const Thyra::TpetraVectorSpace<ST,int,int,NT> > domain = Thyra::tpetraVectorSpace<ST>(Mat->getDomainMap());
   RCP<const Thyra::TpetraVectorSpace<ST,int,int,NT> > range = Thyra::tpetraVectorSpace<ST>(Mat->getRangeMap());

   // convert them to teko compatible sub vectors
   Teko::MultiVector x0_th = Thyra::tpetraVector(domain, x0_tp);
   Teko::MultiVector x1_th = Thyra::tpetraVector(domain, x1_tp);
   Teko::MultiVector b0_th = Thyra::tpetraVector( range, b0_tp);
   Teko::MultiVector b1_th = Thyra::tpetraVector( range, b1_tp);
   std::vector<Teko::MultiVector> x_vec; x_vec.push_back(x0_th); x_vec.push_back(x1_th);
   std::vector<Teko::MultiVector> b_vec; b_vec.push_back(b0_th); b_vec.push_back(b1_th);

   Teko::MultiVector x = Teko::buildBlockedMultiVector(x_vec); // these will be used in the Teko solve
   Teko::MultiVector b = Teko::buildBlockedMultiVector(b_vec);

   // Build the Teko compatible linear system
   Teko::LinearOp thMat = Thyra::tpetraLinearOp<double>(range,domain,Mat);
   Teko::LinearOp thZero;
   Teko::LinearOp A = Thyra::block2x2(thMat,thMat,thZero,thMat); // build an upper triangular 2x2
  

   // Build the preconditioner 
   /////////////////////////////////////////////////////////

   // build an InverseLibrary
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*buildLibPL(),linearSolverBuilder);
   
   // build the inverse factory needed by the example preconditioner
   RCP<Teko::InverseFactory> inverse  = invLib->getInverseFactory("Jacobi");

   // build the preconditioner from the jacobian
   Teko::LinearOp prec = Teko::buildInverse(*inverse,A);

   // Setup the Belos solver
   /////////////////////////////////////////////////////////

   Teuchos::ParameterList belosList;
   belosList.set( "Num Blocks", 200 );               // Maximum number of blocks in Krylov factorization
   belosList.set( "Block Size",1 );                  // Blocksize to be used by iterative solver
   belosList.set( "Maximum Iterations", 200 );       // Maximum number of iterations allowed
   belosList.set( "Maximum Restarts", 1 );           // Maximum number of restarts allowed
   belosList.set( "Convergence Tolerance", 1e-5 );    // Relative convergence tolerance requested
   belosList.set( "Verbosity", 33);//Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
   belosList.set( "Output Frequency", 1 );
   belosList.set( "Output Style", 1 );

   RCP<Belos::LinearProblem<double,MV,OP> > problem = rcp(new Belos::LinearProblem<double,MV,OP>( A, x, b ) );
   problem->setLeftPrec(prec);
   problem->setProblem(); // should check the return type!!!

   RCP<Belos::SolverManager<double,MV,OP> > solver
    = rcp(new Belos::BlockGmresSolMgr<double,MV,OP>(problem, rcpFromRef(belosList)));

   //
   // Perform solve
   //
   Belos::ReturnType ret = solver->solve();

   if (ret!=Belos::Converged) {
     std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
     return -1;
   }

   return 0;
}

RCP<Teuchos::ParameterList> buildLibPL()
{
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   {
      Teuchos::ParameterList & sub = pl->sublist("Jacobi");
      sub.set("Type","Block Jacobi");
      sub.set("Inverse Type","Ifpack2");
   }
   return pl;
}

