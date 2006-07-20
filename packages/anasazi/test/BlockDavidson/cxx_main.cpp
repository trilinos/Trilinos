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
//  This test is for the internal utilities that are used by the modal analysis solver.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBlockDavidson.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicOutputManager.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

#include "ModeLaplace1DQ1.h"

int main(int argc, char *argv[]) 
{
  int i;
  int info = 0;
  
#ifdef EPETRA_MPI
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  Epetra_SerialComm Comm;

#endif
  
  int MyPID = Comm.MyPID();
  
  bool testFailed = false;
  bool verbose = 0;
  std::string which("SM");
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
    else {
      which = argv[1];
    }
  }
  if (argc>2) {
    if (argv[2][0]=='-' && argv[2][1]=='v') {
      verbose = true;
    }
    else {
      which = argv[2];
    }
  }
  
  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }
  
  Anasazi::ReturnType returnCode = Anasazi::Ok;  

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::BasicOutputManager<double>() );

  // Set verbosity level
  if (verbose) {
    MyOM->setVerbosity( Anasazi::FinalSummary + Anasazi::TimingDetails );
  }

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );

  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  // Eigensolver parameters
  int nev = 4;
  int blockSize = 5;
  int maxBlocks = 8;
  int maxIters = 500;
  double tol = 1.0e-6;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Blocks", maxBlocks );
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
  if (info) {
    MyOM->stream(Anasazi::Error) << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;
  }

  // Create the sort manager
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySM = 
     Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );

  // Create the eigensolver  
  Anasazi::BlockDavidson<double, MV, OP> MySolver(MyProblem, MySM, MyOM, MyPL);
  
  // Solve the problem to the specified tolerances or length  
  returnCode = MySolver.solve();
  if (returnCode != Anasazi::Ok)
    testFailed = true;
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();
  Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();
  
  if (verbose)
    info = testCase->eigenCheck( *evecs, &(*evals)[0], 0 );
  
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
  
  for ( i=0; i<nev; i++ ) {
    if ( Teuchos::ScalarTraits<double>::magnitude(normV[i]/(*evals)[i]) > 5.0e-5 )
      testFailed = true;
  }

#ifdef EPETRA_MPI

  MPI_Finalize() ;

#endif

  if (testFailed) {
    if (verbose) {
      MyOM->print(Anasazi::Error,"End Result: TEST FAILED\n");
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose) {
    MyOM->print(Anasazi::Error,"End Result: TEST PASSED\n");
  }
  return 0;
}	
