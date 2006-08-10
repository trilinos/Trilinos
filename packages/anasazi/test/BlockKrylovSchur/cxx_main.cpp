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
// This test is for BlockKrylovSchur solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the BlockKrylovSchurSolMgr solver manager.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"
#include "BlockPCGSolver.h"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{
  bool boolret;
  int MyPID;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;

#endif
  MyPID = Comm.MyPID();

  bool testFailed;
  bool verbose = false;
  bool debug = false;
  std::string filename("mhd1280b.cua");
  std::string which("LM");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  typedef double ScalarType;
  typedef ScalarTraits<ScalarType>                   SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;
  const ScalarType ONE  = SCT::one();

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

  // Create problem
  RefCountPtr<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RefCountPtr<Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RefCountPtr<Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create solver for mass matrix
  const int maxIterCG = 100;
  const double tolCG = 1e-7;
  
  Teuchos::RefCountPtr<BlockPCGSolver> opStiffness = Teuchos::rcp( new BlockPCGSolver(Comm, M.get(), tolCG, maxIterCG, 0) );
  opStiffness->setPreconditioner( 0 );
  Teuchos::RefCountPtr<Anasazi::EpetraGenOp> InverseOp = Teuchos::rcp( new Anasazi::EpetraGenOp( opStiffness, K ) );
  
  // Create the initial vectors
  int blockSize = 5;
  RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create eigenproblem
  const int nev = 5;
  RefCountPtr<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(InverseOp,M,ivec) );
  //
  // Inform the eigenproblem that the operator InverseOp is symmetric under an M inner-product
  problem->setHermitian(true);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  boolret = problem->setProblem();
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
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::StatusTestDetails;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }


  // Eigensolver parameters
  int numBlocks = 10;
  int maxRestarts = 50;
  int stepSize = numBlocks*maxRestarts;
  MagnitudeType tol = tolCG * 10.0;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Step Size", stepSize );
  MyPL.set( "Convergence Tolerance", tol );
  //
  // Create the solver manager
  Anasazi::BlockKrylovSchurSolMgr<ScalarType,MV,OP> MySolverMgr(problem, MyPL);
  // 
  // Check that the parameters were all consumed
  if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
      MyPL.getEntryPtr("Which")->isUsed() == false ||
      MyPL.getEntryPtr("Block Size")->isUsed() == false ||
      MyPL.getEntryPtr("Num Blocks")->isUsed() == false ||
      MyPL.getEntryPtr("Maximum Restarts")->isUsed() == false ||
      MyPL.getEntryPtr("Step Size")->isUsed() == false ||
      MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false) {
    if (verbose && MyPID==0) {
      cout << "Failure! Unused parameters: " << endl;
      MyPL.unused(cout);
    }
  }


  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = problem->getSolution();
  std::vector<MagnitudeType> evals = sol.Evals;
  RefCountPtr<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    ostringstream os;
    os.setf(ios::scientific, ios::floatfield);
    os.precision(6);

    /* finish: this code has bugs: it only works properly when which == "SM" or "SR"
    // Check the problem against the analytical solutions
    if (verbose && which == "LM") {
      info = testCase->eigenCheck( *evecs, &evals[0], 0 );
    }
    */

    // Compute the direct residual
    std::vector<ScalarType> normV( numev );
    SerialDenseMatrix<int,ScalarType> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = evals[i];
    }
    RefCountPtr<MV> Mvecs = MVT::Clone( *evecs, numev ),
                    Kvecs = MVT::Clone( *evecs, numev );
    OPT::Apply( *K, *evecs, *Kvecs );
    OPT::Apply( *M, *evecs, *Mvecs );
    MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
    // compute M-norm of residuals
    OPT::Apply( *M, *Kvecs, *Mvecs );
    MVT::MvDot( *Mvecs, *Kvecs, &normV );
  
    os << "Direct residual norms computed in BlockDavidson_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(M)" << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(evals[i]) != SCT::zero() ) {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) / evals[i] );
      }
      else {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) );
      }
      os << setw(20) << evals[i] << setw(20) << normV[i] << endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (verbose && MyPID==0) {
      cout << endl << os.str() << endl;
    }
  }

#ifdef HAVE_MPI
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
