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
//  Belos as the linear solver within Anasazi.  An Ifpack preconditioner 
//  is constructed to precondition the linear solver.  This operator is 
//  discretized using linear finite elements and constructed as an Epetra 
//  matrix, then passed into the Belos solver to perform the shift-invert
//  operation to be used within the Krylov decomposition.  The specifics 
//  of the block Krylov-Schur method can be set by the user.

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for block Krylov-Schur solver
#include "AnasaziBlockKrylovSchur.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide basic sorting utility required by block Krylov-Schur method
#include "AnasaziBasicSort.hpp"

// Include header to provide Anasazi with Epetra adapters
#include "AnasaziEpetraAdapter.hpp"

// Include header for Epetra compressed-row storage matrix and linear problem
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// Include header for Belos solver and solver interface for Epetra_Operator
#include "BelosEpetraOperator.h"
#include "BelosEpetraAdapter.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestCombo.hpp"

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
  Anasazi::ReturnType returnCode = Anasazi::Ok;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
    Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
  MyOM->SetVerbosity( Anasazi::FinalSummary );  

  std::string which;
  if (argc > 1) {
    which = argv[1];
  }
  else {
    which = "LM";
  }
  if ( which != "SM" && which != "LM" && which != "SR" && which != "LR" ) {
    if (MyOM->doPrint()) {
      std::cout << "Usage: " << argv[0] << " [sort string]" << endl
        << "where:" << endl
        << "sort string       - SM | LM | SR | LR" << endl << endl;
    }
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

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
  if (MyOM->doPrint()) cout << endl << endl;
  if (MyOM->doPrint()) cout << "Constructing ICT preconditioner" << endl;
  int Lfill = 0;
  if (argc > 2) Lfill = atoi(argv[2]);
  if (MyOM->doPrint()) cout << "Using Lfill = " << Lfill << endl;
  int Overlap = 0;
  if (argc > 3) Overlap = atoi(argv[3]);
  if (MyOM->doPrint()) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  if (argc > 4) Athresh = atof(argv[4]);
  if (MyOM->doPrint()) cout << "Using Absolute Threshold Value of " << Athresh << endl;
  double Rthresh = 1.0;
  if (argc >5) Rthresh = atof(argv[5]);
  if (MyOM->doPrint()) cout << "Using Relative Threshold Value of " << Rthresh << endl;
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
  if (MyOM->doPrint()) {
    cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
    cout << endl;
  } 
  //
  //*******************************************************/
  // Set up Belos Block GMRES operator for inner iteration
  //*******************************************************/
  //
  int blockSize = 3;  // block size used by linear solver and eigensolver [ not required to be the same ]
  int maxits = K->NumGlobalRows()/blockSize - 1; // maximum number of iterations to run
  double btol = 1.0e-7;  // relative residual tolerance
  //
  // Create the Belos::LinearProblem
  //
  Teuchos::RefCountPtr<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > 
	  My_LP = Teuchos::rcp( new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>() );
  My_LP->SetOperator( K );
  My_LP->SetLeftPrec( ICT );
  My_LP->SetBlockSize( blockSize );
  //
  // Create the Belos::StatusTest
  //
  Belos::StatusTestMaxIters<double,Epetra_MultiVector,Epetra_Operator> test1( maxits );
  Belos::StatusTestResNorm<double,Epetra_MultiVector,Epetra_Operator> test2( btol );
  Teuchos::RefCountPtr<Belos::StatusTestCombo<double,Epetra_MultiVector,Epetra_Operator> >
	  My_Test = Teuchos::rcp( 
			  new Belos::StatusTestCombo<double,Epetra_MultiVector,Epetra_Operator>
			  ( Belos::StatusTestCombo<double,Epetra_MultiVector,Epetra_Operator>::OR, 
			  test1, test2 ) );
  //
  // Create the Belos::OutputManager
  //
  Teuchos::RefCountPtr<Belos::OutputManager<double> > My_OM = 
	  Teuchos::rcp( new Belos::OutputManager<double>( MyPID ) );
  //My_OM->SetVerbosity( 2 );
  //
  // Create the ParameterList for the Belos Operator
  // 
  Teuchos::RefCountPtr<Teuchos::ParameterList> My_List = Teuchos::rcp( new Teuchos::ParameterList() );
  My_List->set( "Solver", "BlockCG" );
  My_List->set( "MaxIters", maxits );
  //
  // Create the Belos::EpetraOperator
  //
  Teuchos::RefCountPtr<Belos::EpetraOperator> BelosOp = 
    Teuchos::rcp( new Belos::EpetraOperator( My_LP, My_Test, My_OM, My_List ));
  //
  // ************************************
  // Start the block Arnoldi iteration
  // ************************************
  //
  //  Variables used for the Block Arnoldi Method
  //
  int nev = 10;
  int maxBlocks = 3*nev/blockSize;
  int maxRestarts = 5;
  //int step = 5;
  double tol = 1.0e-6;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Max Restarts", maxRestarts );
  MyPL.set( "Tol", tol );
  //MyPL.set( "Step Size", step );
  
  typedef Anasazi::MultiVec<double> MV;
  typedef Anasazi::Operator<double> OP;
  
  // Create an Anasazi::EpetraMultiVec for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> ivec = 
    Teuchos::rcp( new Anasazi::EpetraMultiVec(K->Map(), blockSize) );
  ivec->MvRandom();
  
  // Call the ctor that calls the petra ctor for a matrix
  Teuchos::RefCountPtr<Anasazi::EpetraOp> Kmat = Teuchos::rcp( new Anasazi::EpetraOp(K) );
  Teuchos::RefCountPtr<Anasazi::EpetraOp> Mmat = Teuchos::rcp( new Anasazi::EpetraOp(M) );
  Teuchos::RefCountPtr<Anasazi::EpetraGenOp> Aop = Teuchos::rcp( new Anasazi::EpetraGenOp(BelosOp, M, false) );	
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Aop, Mmat, ivec) );
  
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->SetSymmetric(true);
  
  // Set the number of eigenvalues requested 
  MyProblem->SetNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  info = MyProblem->SetProblem();
  if (info) {
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;
  }

  // Create a sorting manager to handle the sorting of eigenvalues in the solver
  Teuchos::RefCountPtr<Anasazi::BasicSort<double,MV,OP> > MySort = 
    Teuchos::rcp( new Anasazi::BasicSort<double,MV,OP>(which) );
  
  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchur<double,MV,OP> MySolver(MyProblem, MySort, MyOM, MyPL);

  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();

  // Check that the solver returned OK, if not exit example
  if (returnCode != Anasazi::Ok) {
#ifdef EPETRA_MPI
        MPI_Finalize();
#endif
    return -1;
  }
  
  // Obtain eigenvectors directly
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals(); 

  // Retrieve real and imaginary parts of the eigenvectors
  // The size of the eigenvector storage is nev.
  // The real part of the eigenvectors is stored in the first nev vectors.
  // The imaginary part of the eigenvectors is stored in the second nev vectors.
  Anasazi::EpetraMultiVec* evecr = dynamic_cast<Anasazi::EpetraMultiVec*>(MyProblem->GetEvecs()->CloneCopy());

  Teuchos::SerialDenseMatrix<int,double> dmatr(nev,nev);
  Anasazi::EpetraMultiVec tempvec(K->Map(), evecr->GetNumberVecs());	
  K->Apply( *evecr, tempvec );
  tempvec.MvTransMv( 1.0, *evecr, dmatr );

  if (MyOM->doPrint()) {
    double compeval = 0.0;
    cout.setf(ios_base::right, ios_base::adjustfield);
    cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
    cout<<"------------------------------------------------------"<<endl;
    cout<<std::setw(16)<<"Real Part"
	<<std::setw(16)<<"Rayleigh Error"<<endl;
    cout<<"------------------------------------------------------"<<endl;
    for (i=0; i<nev; i++) {
      compeval = dmatr(i,i);
      cout<<std::setw(16)<<compeval
	  <<std::setw(16)<<Teuchos::ScalarTraits<double>::magnitude(compeval-1.0/(*evals)[i])
	  <<endl;
    }
    cout<<"------------------------------------------------------"<<endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}
