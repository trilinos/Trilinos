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
// This example computes the smallest eigenvalues of the discretized 2D Laplacian
// operator using the LOBPCG method.  This problem is discretized using 
// finite elements, resulting in a generalized eigenvalue problem of the form Ax = Mx\lambda.
//
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziLOBPCG.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

#include "ModeLaplace2DQ2.h"

int main(int argc, char *argv[]) 
{
  int i;
  int info = 0;
  
#ifdef EPETRA_MPI
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  
#endif
  
#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  int MyPID = Comm.MyPID();
  
  Anasazi::ReturnType returnCode = Anasazi::Ok;  

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  // Number of dimension of the domain
  int space_dim = 2;

  // Size of each of the dimensions of the domain
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  brick_dim[1] = 1.0;

  // Number of elements in each of the dimensions of the domain
  std::vector<int> elements( space_dim );
  elements[0] = 10;
  elements[1] = 10;
  
  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );

  // Set verbosity level
  MyOM->SetVerbosity( Anasazi::FinalSummary );

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );
  
  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_Operator> K = Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_Operator> M = Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getMass()), false );
  
  // Eigensolver parameters
  int nev = 10;
  int blockSize = 5;
  int maxIters = 500;
  double tol = 1.0e-6;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Iters", maxIters );
  MyPL.set( "Tol", tol );
  
  // Create eigenproblem

  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(K, M, ivec) );
  
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true);
  
  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->SetNEV( nev );
  
  // Inform the eigenproblem that you are finishing passing it information
  info = MyProblem->SetProblem();
  if (info)
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;
  
  // Create the eigensolver
  
  Anasazi::LOBPCG<double, MV, OP> MySolver(MyProblem, MyOM, MyPL);
  
  // Solve the problem to the specified tolerances or length

  returnCode = MySolver.solve();

  // Check that the solver returned Ok, if not exit example
  if (returnCode != Anasazi::Ok)
    return -1;
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();
  Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();
  
  // Compute the direct residual
  std::vector<double> normV( evecs->NumVectors() );
  Teuchos::SerialDenseMatrix<int,double> T(evecs->NumVectors(), evecs->NumVectors());
  for (int i=0; i<evecs->NumVectors(); i++)
    T(i,i) = (*evals)[i];
  Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
  K->Apply( *evecs, Kvec );  
  Epetra_MultiVector Mvec( M->OperatorDomainMap(), evecs->NumVectors() );
  M->Apply( *evecs, Mvec );  
  Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
  info = Kvec.Norm2( &normV[0] );
  assert( info==0 );

  if (MyPID == 0) {
    cout<<"Actual Residuals"<<endl;
    cout<<"------------------------------------------------------"<<endl;
    cout<<"Eigenvalue"<<"\t\t"<<"Direct Residual"<<endl;
    cout<<"------------------------------------------------------"<<endl;
    for (i=0; i<nev; i++) {
      cout<< (*evals)[i] << "\t\t"<< normV[i]/(*evals)[i] << endl;
    }
    cout<<"------------------------------------------------------"<<endl;
  }
  //
  // Default return value
  //
  return 0;
}
