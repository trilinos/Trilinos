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

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

#include "ModeLaplace1DQ1.h"
#include "BlockPCGSolver.h"

int main(int argc, char *argv[]) 
{
  int i, info = 0;
  double zero = 0.0;
  
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
  int NumProc = Comm.NumProc();
  
  bool testFailed = false;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  
  if (verbose && MyPID == 0)
    cout << Anasazi::Anasazi_Version() << endl << endl;
  
  int numberFailedTests = 0;
  int returnCode = 0;  

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 50;

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::OutputManager<double>() );

  // Set verbosity level
  if (verbose && MyPID == 0)
    MyOM->SetVerbosity( Anasazi::FinalSummary );

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );

  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_Operator> K = Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_Operator> M = Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getMass()), false );

  // Create preconditioner
  int maxIterCG = 100;
  double tolCG = 1e-05;
  
  Teuchos::RefCountPtr<BlockPCGSolver> opStiffness = Teuchos::rcp( new BlockPCGSolver(Comm, K.get(), tolCG, maxIterCG, 3) );
  opStiffness->setPreconditioner( 0 );

  int nev = 4;
  int blockSize = 5;
  int maxBlocks = 8;
  int maxIter = 500;
  double tol = tolCG * 10.0;
  
  // Create eigenproblem

  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(opStiffness, M, ivec) );
  MyProblem->SetPrec( Teuchos::rcp( const_cast<Epetra_Operator *>(opStiffness->getPreconditioner()), false ) );
  
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true);
  
  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->SetNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  assert( MyProblem->SetProblem() == 0 );

  // Create the eigensolver

  Anasazi::BlockDavidson<double, MV, OP> MySolver(MyProblem, MyOM, tol,
						  blockSize, maxBlocks, maxIter);

  // Solve the problem to the specified tolerances or length
  MySolver.solve();

  return 0;
}	
