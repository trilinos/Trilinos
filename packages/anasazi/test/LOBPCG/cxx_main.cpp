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
// This test is for LOBPCG solving a standard (Ax=xl) real Hermitian
// eigenvalue problem, using the LOBPCGSolMgr solver manager.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"
using namespace Teuchos;

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

  bool testFailed;
  bool verbose = 0;
  std::string which("LM");
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

  Anasazi::ReturnType returnCode;
  
  typedef double ScalarType;
  typedef ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

  RefCountPtr<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RefCountPtr<Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RefCountPtr<Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create the initial vectors
  int blockSize = 5;
  RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create eigenproblem
  int nev = 5;
  RefCountPtr<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > MyProblem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(K,M,ivec) );
  //
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->setHermitian(true);
  //
  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
           << "End Result: TEST FAILED" << endl;	
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }


  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }


  // Eigensolver parameters
  int maxIters = 450;
  MagnitudeType tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10 );
  MyPL.set( "Which", which );
  MyPL.set( "Full Ortho", true );
  MyPL.set( "Verbosity", verbosity );
  //
  // Create the solver manager
  Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  returnCode = MySolverMan.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  std::vector<MagnitudeType> evals = sol.Evals;
  RefCountPtr<Epetra_MultiVector> evecs = sol.Evecs;
  
  // Check the problem against the analytical solutions
  if (verbose && returnCode == Anasazi::Converged) {
    info = testCase->eigenCheck( *evecs, &evals[0], 0 );
  }
  
  // Compute the direct residual
  std::vector<MagnitudeType> normV( evecs->NumVectors() );
  SerialDenseMatrix<int,ScalarType> T(evecs->NumVectors(), evecs->NumVectors());
  for (int i=0; i<evecs->NumVectors(); i++) {
    T(i,i) = evals[i];
  }
  Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
  K->Apply( *evecs, Kvec );  
  Epetra_MultiVector Mvec( M->OperatorDomainMap(), evecs->NumVectors() );
  M->Apply( *evecs, Mvec );  
  Anasazi::MultiVecTraits<ScalarType,Epetra_MultiVector>::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
  info = Kvec.Norm2( &normV[0] );
  assert( info==0 );
  
  for ( i=0; i<nev; i++ ) {
    if ( ScalarTraits<ScalarType>::magnitude(normV[i]/evals[i]) > 5.0e-5 ) {
      testFailed = true;
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}	
