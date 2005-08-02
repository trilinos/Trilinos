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
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "ModeLaplace1DQ1.h"
#include "BlockPCGSolver.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

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
  
  Anasazi::ReturnType returnCode = Anasazi::Ok;

  bool testFailed = false;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  
  if (verbose && MyPID == 0)
    cout << Anasazi::Anasazi_Version() << endl << endl;
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );

  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_Operator> K = Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_Operator> M = Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getMass()), false );

  // Create preconditioner
  int maxIterCG = 100;
  double tolCG = 1e-7;
  
  Teuchos::RefCountPtr<BlockPCGSolver> opStiffness = Teuchos::rcp( new BlockPCGSolver(Comm, M.get(), tolCG, maxIterCG, 0) );
  opStiffness->setPreconditioner( 0 );

  Teuchos::RefCountPtr<Anasazi::EpetraGenOp> InverseOp = Teuchos::rcp( new Anasazi::EpetraGenOp( opStiffness, K ) );

  // Variables needed for eigensolver

  int nev = 4;
  int blockSize = 2;
  int maxBlocks = 10;
  int maxRestarts = 500;
  double tol = tolCG * 10.0;
  std::string which = "SM";  
  
  // Create parameter list to pass into solver

  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Max Restarts", maxRestarts );
  MyPL.set( "Tol", tol );
  
  // Create eigenproblem
  
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(InverseOp, M, ivec) );
  
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true);
  
  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->SetNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  info = MyProblem->SetProblem();
  if (info)
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );

  // Set verbosity level
  if (verbose) 
    MyOM->SetVerbosity( Anasazi::FinalSummary + Anasazi::TimingDetails );

  // Create a sorting manager to handle the sorting of eigenvalues in the solver
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
    Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );
  
  // Create the eigensolver  
  Anasazi::BlockKrylovSchur<double, MV, OP> MySolver(MyProblem, MySort, MyOM, MyPL);
  
  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();
  if (returnCode != Anasazi::Ok) {
    cout << "The return code for BlockKrylovSchur is : "<< returnCode << endl;
    testFailed = true;
  }
  
  // Obtain results directly
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();
  Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();
  
  // Check eigenvalues with ModeLaplace
  if (verbose)
    info = testCase->eigenCheck( *evecs, &(*evals)[0], 0 );

  Epetra_MultiVector tempvec(K->OperatorDomainMap(), evecs->NumVectors());	
  K->Apply( *evecs, tempvec );
  Epetra_MultiVector Mtempvec(M->OperatorDomainMap(), evecs->NumVectors());	
  M->Apply( *evecs, Mtempvec );

  Epetra_LocalMap LocalMap(nev, 0, Comm);
  Epetra_MultiVector dmatr( LocalMap, nev );

  for ( i=0; i<nev; i++ ) {
    dmatr[i][i] = (*evals)[i];
  }
  tempvec.Multiply( 'N', 'N', -1.0, Mtempvec, dmatr, 1.0 );
  std::vector< double > normvec( nev );
  tempvec.Norm2( &normvec[0] );
  
  for ( i=0; i<nev; i++ ) {
    if ( Teuchos::ScalarTraits<double>::magnitude(normvec[i]/(*evals)[i]) > 5.0e-5 )
      testFailed = true;
  }

#ifdef EPETRA_MPI

  MPI_Finalize() ;

#endif

  if (testFailed) {
    if (verbose && MyPID==0)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;

}	
