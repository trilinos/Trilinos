// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
//  This example computes the specified eigenvalues of the discretized 2D Laplacian
//  using the block Davidson.  

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockJacobiDavidson.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"

typedef double ScalarType;
typedef double MagnitudeType;

#include "MyOperator.hpp"
#include "MyMultiVec.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) 
{
  int i, info;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  Anasazi::ReturnType returnCode = Anasazi::Ok;	

  int            DIM  = 100;
  ScalarType     TOL    = 1e-8;

  int MAXITER = 500;
  int NEV = 3; 
  int SMIN = 20;
  int SMAX = 30;
  // FIXME: if I set BLOCKSIZE = 2 then nothing good happens...
  int BLOCKSIZE = 2;

  // ============================================= //
  // Sets up the test problem, using the operators //
  // A = tridiag([-1,4,-1]), B=tridiag([-1,2,-1])  //
  // ============================================= //

  ScalarType TARGET = 1.5;

  std::vector<ScalarType> A_Entries(3);
  A_Entries[0] = -1.0;
  A_Entries[1] =  4.0;
  A_Entries[2] = -1.0;
  Teuchos::RefCountPtr<MyOperator> A = Teuchos::rcp(new MyOperator(DIM, A_Entries));

  std::vector<ScalarType> B_Entries(3);
  B_Entries[0] = -1.0;
  B_Entries[1] =  2.0;
  B_Entries[2] = -1.0;
  Teuchos::RefCountPtr<MyOperator> B = Teuchos::rcp(new MyOperator(DIM, B_Entries));

  std::vector<ScalarType> K_Entries(3);
  K_Entries[0] =  0.0;
  K_Entries[1] =  1.0 / (A_Entries[1] - TARGET * B_Entries[1]);
  K_Entries[2] =  0.0;
  Teuchos::RefCountPtr<MyOperator> K = Teuchos::rcp(new MyOperator(DIM, K_Entries));
  
  // ================== //
  // Sets up the solver //
  // ================== //

  typedef Anasazi::MultiVec<ScalarType> MV;        
  typedef Anasazi::Operator<ScalarType> OP;        
  //typedef Anasazi::MultiVecTraits<ScalarType, MyMultiVec> MVT;

  Teuchos::RefCountPtr<MV> ivec = Teuchos::rcp(new MyMultiVec(DIM, BLOCKSIZE));        
  ivec->MvRandom();

  // Create the eigenproblem.
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<ScalarType, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<ScalarType, MV, OP>(A, B, ivec) );

  MyProblem->SetPrec(K); 

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true); 

  // Set the number of eigenvalues requested
  MyProblem->SetNEV(NEV);

  // Inform the eigenproblem that you are finishing passing it information
  MyProblem->SetProblem();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<ScalarType> > MyOM =
    Teuchos::rcp(new Anasazi::OutputManager<ScalarType>(MyPID));
  MyOM->SetVerbosity(Anasazi::FinalSummary);	

  // Create a sort manager
  Teuchos::RefCountPtr<Anasazi::BasicSort<ScalarType, MV, OP> > MySM =
    Teuchos::rcp(new Anasazi::BasicSort<ScalarType, MV, OP>("SM"));


  // Create parameter list to pass into solver
  // FIXME: ADD PARAMTERS
  Teuchos::ParameterList MyPL;
  MyPL.set("Block Size", BLOCKSIZE);
  MyPL.set("SMIN", SMIN);
  MyPL.set("SMAX", SMAX);
  MyPL.set("Max Iters", MAXITER);
  MyPL.set("Tol", TOL);
  MyPL.set("Target", TARGET);

  // Initialize the Block Jacobi-Davidson solver
  Anasazi::BlockJacobiDavidson<ScalarType, MagnitudeType, MV, OP> MySolver(MyProblem, MySM, MyOM, MyPL);
                           
  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();

  MySolver.currentStatus();

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  //
  // Default return value
  //
  return(EXIT_SUCCESS);
}
