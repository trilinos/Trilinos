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
//  This example computes the eigenvalues of smallest magnitude of the 
//  discretized 2D Laplacian operator using the block Krylov-Schur method.  
//  This problem shows the construction of an inner-outer iteration using 
//  AztecOO as the linear solver within Anasazi.  An Ifpack preconditioner 
//  is constructed to precondition the linear solver.  This operator is 
//  discretized using linear finite elements and constructed as an Epetra 
//  matrix, then passed into the AztecOO solver to perform the shift-invert
//  operation to be used within the Krylov decomposition.  The specifics 
//  of the block Krylov-Schur method can be set by the user.

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for block Krylov-Schur solver manager
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide Anasazi with Epetra adapters
#include "AnasaziEpetraAdapter.hpp"

// Include header for Epetra compressed-row storage matrix and linear problem
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// Include header for AztecOO solver and solver interface for Epetra_Operator
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// Include header for Ifpack incomplete Cholesky preconditioner
#include "Ifpack_CrsIct.h"

// Include header for Teuchos serial dense matrix
#include "Teuchos_SerialDenseMatrix.hpp"

// Include header for the problem definition
#include "ModeLaplace2DQ2.h"

// Include selected communicator class and map required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  int i, info;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

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
  
  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );
  
  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  
  //
  //*****Select the Preconditioner*****
  //
  if (MyPID==0) cout << endl << endl;
  if (MyPID==0) cout << "Constructing ICT preconditioner" << endl;
  int Lfill = 0;
  if (argc > 2) Lfill = atoi(argv[2]);
  if (MyPID==0) cout << "Using Lfill = " << Lfill << endl;
  int Overlap = 0;
  if (argc > 3) Overlap = atoi(argv[3]);
  if (MyPID==0) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  if (argc > 4) Athresh = atof(argv[4]);
  if (MyPID==0) cout << "Using Absolute Threshold Value of " << Athresh << endl;
  double Rthresh = 1.0;
  if (argc >5) Rthresh = atof(argv[5]);
  if (MyPID==0) cout << "Using Relative Threshold Value of " << Rthresh << endl;
  double dropTol = 1.0e-6;
  //
  Teuchos::RefCountPtr<Ifpack_CrsIct> ICT;
  //
  if (Lfill > -1) {
    ICT = Teuchos::rcp( new Ifpack_CrsIct(*K, dropTol, Lfill) );
    ICT->SetAbsoluteThreshold(Athresh);
    ICT->SetRelativeThreshold(Rthresh);
    int initerr = ICT->InitValues(*K);
    if (initerr != 0) cout << "InitValues error = " << initerr;
    info = ICT->Factor();
    assert( info==0 );
  } 
  //
  bool transA = false;
  double Cond_Est;
  ICT->Condest(transA, Cond_Est);
  if (MyPID==0) {
    cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
    cout << endl;
  } 
  //
  // *******************************************************
  // Set up AztecOO GMRES operator for inner iteration
  // *******************************************************
  //
  // Create Epetra linear problem class to solve "Ax = b"
  Epetra_LinearProblem precProblem;
  precProblem.SetOperator(K.get());
  
  // Create AztecOO solver for solving "Ax = b" using an incomplete cholesky preconditioner
  AztecOO precSolver(precProblem);
  precSolver.SetPrecOperator(ICT.get());
  precSolver.SetAztecOption(AZ_output, AZ_none);
  precSolver.SetAztecOption(AZ_solver, AZ_cg);
  
  // Use AztecOO solver to create the AztecOO_Operator
  Teuchos::RefCountPtr<AztecOO_Operator> precOperator =
    Teuchos::rcp( new AztecOO_Operator(&precSolver, 100) );	
  //
  // ************************************
  // Start the block Arnoldi iteration
  // ************************************
  //
  //  Variables used for the Block Arnoldi Method
  //
  int nev = 10;
  int blockSize = 3;  
  int numBlocks = 3*nev/blockSize;
  int maxRestarts = 5;
  //int step = 5;
  double tol = 1.0e-8;
  std::string which = "LM";
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  //MyPL.set( "Step Size", step );
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
  
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->Map(), blockSize) );
  MVT::MvRandom( *ivec );
  
  // Call the ctor that calls the petra ctor for a matrix
  Teuchos::RefCountPtr<Anasazi::EpetraGenOp> Aop = Teuchos::rcp( new Anasazi::EpetraGenOp(precOperator, M) );	
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Aop, M, ivec) );
  
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian(true);
  
  // Set the number of eigenvalues requested 
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RefCountPtr<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    Teuchos::RefCountPtr<MV> evecs = sol.Evecs;
    int numev = sol.numVecs;
    
    Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
    Epetra_MultiVector tempvec(K->Map(), MVT::GetNumberVecs( *evecs ));	
    OPT::Apply( *K, *evecs, tempvec );
    MVT::MvTransMv( 1.0, tempvec, *evecs, dmatr );
    
    if (MyPID==0) {
      double compeval = 0.0;
      cout.setf(ios_base::right, ios_base::adjustfield);
      cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
      cout<<"------------------------------------------------------"<<endl;
      cout<<std::setw(16)<<"Real Part"
	  <<std::setw(16)<<"Rayleigh Error"<<endl;
      cout<<"------------------------------------------------------"<<endl;
      for (i=0; i<numev; i++) {
	compeval = dmatr(i,i);
	cout<<std::setw(16)<<compeval
	    <<std::setw(16)<<Teuchos::ScalarTraits<double>::magnitude(compeval-1.0/evals[i].realpart)
	    <<endl;
      }
      cout<<"------------------------------------------------------"<<endl;
    }
    
  }
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}
