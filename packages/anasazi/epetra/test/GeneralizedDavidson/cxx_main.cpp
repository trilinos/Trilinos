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
// This test is for GeneralizedDavidson solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the GeneralizedDavidsonSolMgr solver manager.
// A few steps of BlockPCG is used as a preconditioner.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_InvOperator.h"

#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
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
  if (debug) verbose = true;

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
  RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RCP<Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RCP<Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create solver for mass matrix
  // Note that accuracy of Davidson solution does NOT depend on how accurately the BlockPCG is solved
  const int maxIterCG = 10;
  const double tolCG = 1e-2;

  RCP<BlockPCGSolver> opStiffness = rcp( new BlockPCGSolver(Comm, M.get(), tolCG, maxIterCG, 0) );
  opStiffness->setPreconditioner( 0 );
  RCP<Epetra_Operator> invStiffness = rcp( new Epetra_InvOperator(opStiffness.get()) );

  // Create the initial vectors
  int blockSize = 3;
  RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create eigenproblem
  const int nev = 5;
  RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>() );
  problem->setA(K);
  problem->setM(M);
  problem->setPrec(invStiffness);
  problem->setInitVec(ivec);
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
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }


  // Eigensolver parameters
  int maxRestarts = 25;
  int maxDim = 50;
  MagnitudeType tol = 1e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Maximum Subspace Dimension", maxDim );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  //
  // Create the solver manager
  Anasazi::GeneralizedDavidsonSolMgr<ScalarType,MV,OP> MySolverMgr(problem, MyPL);
  //
  // Check that the parameters were all consumed
  if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
      MyPL.getEntryPtr("Which")->isUsed() == false ||
      MyPL.getEntryPtr("Maximum Subspace Dimension")->isUsed() == false ||
      MyPL.getEntryPtr("Block Size")->isUsed() == false ||
      MyPL.getEntryPtr("Maximum Restarts")->isUsed() == false ||
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
  std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<ScalarType> normV( numev );
    SerialDenseMatrix<int,ScalarType> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = evals[i].realpart;
    }
    RCP<MV> Mvecs = MVT::Clone( *evecs, numev );
    RCP<MV> Kvecs = MVT::Clone( *evecs, numev );
    OPT::Apply( *K, *evecs, *Kvecs );
    OPT::Apply( *M, *evecs, *Mvecs );
    MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
    // compute 2-norm of residuals
    std::vector<MagnitudeType> resnorm(numev);
    MVT::MvNorm( *Kvecs, resnorm );

    os << "Number of iterations performed in GeneralizedDavidson_test.exe: " << MySolverMgr.getNumIters() << endl
       << "Direct residual norms computed in GeneralizedDavidson_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual" << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      os << std::setw(20) << evals[i].realpart << std::setw(20) << resnorm[i] << endl;
      if ( resnorm[i] > tol ) {
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
