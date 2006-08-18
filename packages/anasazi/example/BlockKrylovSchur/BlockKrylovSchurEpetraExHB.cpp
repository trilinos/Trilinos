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
//  This example compute the eigenvalues of a Harwell-Boeing matrix using the block Arnoldi
//  method.  The matrix is passed to the example routine through the command line, and 
//  converted to an Epetra matrix through some utilty routines.  This matrix is passed to the
//  eigensolver and then used to construct the Krylov decomposition.  The specifics of the 
//  block Arnoldi method can be set by the user.

#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "Trilinos_Util.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  //
  int i,info;
  int n_nonzeros, N_update;
  int *bindx=0, *update=0;
  double *val=0;

#ifdef EPETRA_MPI  
  // Initialize MPI  
  MPI_Init(&argc,&argv);   
  Epetra_MpiComm Comm( MPI_COMM_WORLD );  
#else  
  Epetra_SerialComm Comm;  
#endif
  
  int MyPID = Comm.MyPID();

  bool verbose=false;
  std::string filename = "";
  std::string which = "LM";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,or LR).");
  cmdp.setOption("filename",&filename,"Filename and path of a Harwell-Boeing data set.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  //
  //**********************************************************************
  //******************Set up the problem to be solved*********************
  //**********************************************************************
  //
  int NumGlobalElements;  // total # of rows in matrix
  //
  // *****Read in matrix from HB file******
  //
  Trilinos_Util_read_hb(const_cast<char *>(filename.c_str()), MyPID, &NumGlobalElements, 
			&n_nonzeros, &val, &bindx);
  //
  // *****Distribute data among processors*****
  //
  Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalElements, &n_nonzeros, &N_update,
                                 &update, &val, &bindx);
  //
  // *****Construct the matrix*****
  //
  int NumMyElements = N_update; // # local rows of matrix on processor
  //
  // Create an integer vector NumNz that is used to build the Epetra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  //
  std::vector<int> NumNz(NumMyElements);
  for (i=0; i<NumMyElements; i++) {
    NumNz[i] = bindx[i+1] - bindx[i] + 1;
  }
  //
  Epetra_Map Map(NumGlobalElements, NumMyElements, update, 0, Comm);
  //
  // Create a Epetra_Matrix
  //
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]) );
  //
  // Add rows one-at-a-time
  //
  int NumEntries;
  for (i=0; i<NumMyElements; i++) {
    NumEntries = bindx[i+1] - bindx[i];
    info = A->InsertGlobalValues(update[i], NumEntries, val + bindx[i], bindx + bindx[i]);
    assert(info==0 );
    info = A->InsertGlobalValues(update[i], 1, val+i, update+i);
    assert( info==0 );
  }
  //
  // Finish up
  //
  info = A->FillComplete();
  assert( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  //
  //************************************
  // Start the block Arnoldi iteration
  //***********************************         
  //
  //  Variables used for the Block Arnoldi Method
  // 
  int nev = 5;
  int blockSize = 5;
  int numBlocks = 10;
  int maxRestarts = 10;
  int step = 5;
  double tol = 1.0e-8;
  
  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
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
  MyPL.set( "Step Size", step );

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
  //
  // Create the eigenproblem to be solved.
  //
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );
  
  // Inform the eigenproblem that the matrix A is Hermitian
  //MyProblem->setHermitian(true);

  // Set the number of eigenvalues requested 
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
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
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RefCountPtr<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;

  if (numev > 0) {
    // Compute residuals.
    Teuchos::LAPACK<int,double> lapack;
    std::vector<double> normA(numev);
    
    if (MyProblem->isHermitian()) {
      // Get storage
      Epetra_MultiVector Aevecs(Map,numev);
      Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
      B.putScalar(0.0); 
      for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}
      
      // Compute A*evecs
      OPT::Apply( *A, *evecs, Aevecs );
      
      // Compute A*evecs - lambda*evecs and its norm
      MVT::MvTimesMatAddMv( -1.0, *evecs, B, 1.0, Aevecs );
      MVT::MvNorm( Aevecs, &normA );
      
      // Scale the norms by the eigenvalue
      for (int i=0; i<numev; i++) {
	normA[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
      }
    } else {
      // The problem is non-Hermitian.
      int i=0;
      std::vector<int> curind(1);
      std::vector<double> resnorm(1), tempnrm(1);
      Teuchos::RefCountPtr<MV> evecr, eveci, tempAevec;
      Epetra_MultiVector Aevec(Map,numev);
      
      // Compute A*evecs
      OPT::Apply( *A, *evecs, Aevec );
      
      Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
      while (i<numev) {
	if (index[i]==0) {
	  // Get a view of the current eigenvector (evecr)
	  curind[0] = i;
	  evecr = MVT::CloneView( *evecs, curind );
	  
	  // Get a copy of A*evecr
	  tempAevec = MVT::CloneCopy( Aevec, curind );
	  
	  // Compute A*evecr - lambda*evecr
	  Breal(0,0) = evals[i].realpart;
	  MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
	  
	  // Compute the norm of the residual and increment counter
	  MVT::MvNorm( *tempAevec, &resnorm );
	  normA[i] = resnorm[0]/Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
	  i++;
	} else {
	  // Get a view of the real part of the eigenvector (evecr)
	  curind[0] = i;
	  evecr = MVT::CloneView( *evecs, curind );
	  
	  // Get a copy of A*evecr
	  tempAevec = MVT::CloneCopy( Aevec, curind );
	  
	  // Get a view of the imaginary part of the eigenvector (eveci)
	  curind[0] = i+1;
	  eveci = MVT::CloneView( *evecs, curind );
	  
	  // Set the eigenvalue into Breal and Bimag
	  Breal(0,0) = evals[i].realpart;
	  Bimag(0,0) = evals[i].imagpart;
	  
	  // Compute A*evecr - evecr*lambdar + eveci*lambdai
	  MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
	  MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
	  MVT::MvNorm( *tempAevec, &tempnrm );
	  
	  // Get a copy of A*eveci
	  tempAevec = MVT::CloneCopy( Aevec, curind );
	  
	  // Compute A*eveci - eveci*lambdar - evecr*lambdai
	  MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
	  MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
	  MVT::MvNorm( *tempAevec, &resnorm );
	  
	  // Compute the norms and scale by magnitude of eigenvalue
	  normA[i] = lapack.LAPY2( tempnrm[i], resnorm[i] ) /
	    lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
	  normA[i+1] = normA[i];
	  
	  i=i+2;
	}
      }
    }
    
    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      cout.setf(ios_base::right, ios_base::adjustfield);	
      cout<<endl<< "Actual Residuals"<<endl;
      if (MyProblem->isHermitian()) {
	cout<< std::setw(16) << "Real Part"
	    << std::setw(20) << "Direct Residual"<< endl;
	cout<<"-----------------------------------------------------------"<<endl;
	for (int i=0; i<numev; i++) {
	  cout<< std::setw(16) << evals[i].realpart 
	      << std::setw(20) << normA[i] << endl;
	}  
	cout<<"-----------------------------------------------------------"<<endl;
      } 
      else {
	cout<< std::setw(16) << "Real Part"
	    << std::setw(16) << "Imag Part"
	    << std::setw(20) << "Direct Residual"<< endl;
	cout<<"-----------------------------------------------------------"<<endl;
	for (int i=0; i<numev; i++) {
	  cout<< std::setw(16) << evals[i].realpart 
	      << std::setw(16) << evals[i].imagpart 
	      << std::setw(20) << normA[i] << endl;
	}  
	cout<<"-----------------------------------------------------------"<<endl;
      }  
    }
  }

  if (bindx) delete [] bindx;
  if (update) delete [] update;
  if (val) delete [] val;

#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return 0;

} // end BlockKrylovSchurEpetraExHb.cpp
