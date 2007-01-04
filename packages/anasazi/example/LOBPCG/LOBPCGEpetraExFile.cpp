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
//  This example compute the eigenvalues of a matrix from and input file using the block Davidson 
//  method.  The matrix is passed to the example routine through the command line, and 
//  converted to an Epetra matrix through some utilty routines.  This matrix is passed to the
//  eigensolver, the specifics of the block Davidson method can be set by the user.

#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  //
  bool haveM = false;

#ifdef EPETRA_MPI  
  // Initialize MPI  
  MPI_Init(&argc,&argv);   
  Epetra_MpiComm Comm( MPI_COMM_WORLD );  
#else  
  Epetra_SerialComm Comm;  
#endif
  
  int MyPID = Comm.MyPID();

  bool verbose=false;
  std::string k_filename = "";
  std::string m_filename = "";
  std::string which = "SM";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,or LR).");
  cmdp.setOption("K-filename",&k_filename,"Filename and path of the stiffness matrix.");
  cmdp.setOption("M-filename",&m_filename,"Filename and path of the mass matrix.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (k_filename=="") {
    cout << "The matrix K must be supplied through an input file!!!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (m_filename!="") {
    haveM = true;
  }
  //
  //**********************************************************************
  //******************Set up the problem to be solved*********************
  //**********************************************************************
  //
  // *****Read in matrix from file******
  //
  Teuchos::RefCountPtr<Epetra_Map> Map;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> K, M;
  EpetraExt::readEpetraLinearSystem( k_filename, Comm, &K, &Map );

  if (haveM) {
    EpetraExt::readEpetraLinearSystem( m_filename, Comm, &M, &Map );
  }
  //
  //************************************
  // Start the block Davidson iteration 
  //***********************************         
  //
  //  Variables used for the LOBPCG Method
  // 
  int nev = 5;
  int blockSize = 5;
  int maxRestarts = 1000;
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
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
  //
  // Create the eigenproblem to be solved.
  //
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(*Map, blockSize) );
  ivec->Random();
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem;
  if (haveM) {
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>( K, M, ivec ) );
  } 
  else {
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>( K, ivec ) );
  } 
 
  // Inform the eigenproblem that (K,M) is Hermitian
  MyProblem->setHermitian(true);

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

  // Initialize the LOBPCG solver
  Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);
    
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
    std::vector<double> normEV(numev);
    
    // Get storage
    Teuchos::RefCountPtr<Epetra_MultiVector> Kevecs, Mevecs;
    Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
    B.putScalar(0.0); 
    for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}
    
    // Compute K*evecs
    Kevecs = Teuchos::rcp(new Epetra_MultiVector(*Map,numev) );
    OPT::Apply( *K, *evecs, *Kevecs );

    // Compute M*evecs
    if (haveM) {
      Mevecs = Teuchos::rcp(new Epetra_MultiVector(*Map,numev) );
      OPT::Apply( *M, *evecs, *Mevecs );
    }
    else {
      Mevecs = evecs;
    }

    // Compute K*evecs - lambda*M*evecs and its norm
    MVT::MvTimesMatAddMv( -1.0, *Mevecs, B, 1.0, *Kevecs );
    MVT::MvNorm( *Kevecs, &normEV );
    
    // Scale the norms by the eigenvalue
    for (int i=0; i<numev; i++) {
      normEV[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
    }
    
    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      cout.setf(ios_base::right, ios_base::adjustfield);	
      cout<<endl<< "Actual Residuals"<<endl;
      cout<< std::setw(16) << "Real Part"
	  << std::setw(20) << "Direct Residual"<< endl;
      cout<<"-----------------------------------------------------------"<<endl;
      for (int i=0; i<numev; i++) {
	cout<< std::setw(16) << evals[i].realpart 
	    << std::setw(20) << normEV[i] << endl;
      }  
      cout<<"-----------------------------------------------------------"<<endl;
    }
  }
  
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  return 0;
  
} // end LOBPCGEpetraExFile.cpp
