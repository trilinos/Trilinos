// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This test is for LOBPCG solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the LOBPCGSolMgr solver manager.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziSVQBOrthoManager.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

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
  bool shortrun = false;
  bool fullOrtho = true;
  bool testRecovery = false;
  std::string which("LM");
  int numElements = 100;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
  cmdp.setOption("testrecovery","notestrecovery",&testRecovery,"Test the LOBPCGRitzError recovery code in LOBPCGSolMgr.");
  cmdp.setOption("fullortho","nofullortho",&fullOrtho,"Use full orthogonalization.");
  cmdp.setOption("numElements",&numElements,"Number of elements in discretization.");
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
    std::cout << Anasazi::Anasazi_Version() << std::endl << std::endl;
  }

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = numElements;

  // Create problem
  RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RCP<const Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create the initial vectors
  int blockSize = 5;
  RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create eigenproblem
  const int nev = 5;
  RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(K,M,ivec) );
  //
  // Inform the eigenproblem that the operator K is symmetric
  problem->setHermitian(true);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  boolret = problem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      std::cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << std::endl
           << "End Result: TEST FAILED" << std::endl;	
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
  int maxIters;
  if (shortrun) {
    maxIters = 50;
  }
  else {
    maxIters = 450;
  }
  MagnitudeType tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10 );
  MyPL.set( "Full Ortho", fullOrtho );
  if (testRecovery) {
    // initialize with bad P
    RCP<Anasazi::LOBPCGState<ScalarType,MV> > badstate = rcp( new Anasazi::LOBPCGState<ScalarType,MV>() );
    Anasazi::SVQBOrthoManager<ScalarType,MV,OP> ortho(M);
    ortho.normalize(*ivec,null);
    badstate->X = ivec;
    RCP<MV> P  = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
    MVT::MvInit(*P,0.0);
    badstate->P = P;
    MyPL.set( "Init", badstate );
  }
  //
  // Create the solver manager
  Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> MySolverMan(problem, MyPL);
  // 
  // Check that the parameters were all consumed
  if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
      MyPL.getEntryPtr("Which")->isUsed() == false ||
      MyPL.getEntryPtr("Block Size")->isUsed() == false ||
      MyPL.getEntryPtr("Maximum Iterations")->isUsed() == false ||
      MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false ||
      MyPL.getEntryPtr("Use Locking")->isUsed() == false ||
      MyPL.getEntryPtr("Locking Tolerance")->isUsed() == false ||
      MyPL.getEntryPtr("Full Ortho")->isUsed() == false) {
    if (verbose && MyPID==0) {
      std::cout << "Failure! Unused parameters: " << std::endl;
      MyPL.unused(std::cout);
    }
  }


  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
    if ( shortrun==false && (testRecovery==false || fullOrtho==false) ) {
      testFailed = true;
    }
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

    // Check the problem against the analytical solutions
    if (verbose) {
      double *revals = new double[numev];
      for (int i=0; i<numev; i++) {
        revals[i] = evals[i].realpart;
      }
      bool smallest = false;
      if (which == "SM" || which == "SR") {
        smallest = true;
      }
      testCase->eigenCheck( *evecs, revals, 0, smallest );
      delete [] revals;
    }

    // Compute the direct residual
    std::vector<ScalarType> normV( numev );
    SerialDenseMatrix<int,ScalarType> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = evals[i].realpart;
    }
    RCP<MV> Mvecs = MVT::Clone( *evecs, numev ),
                    Kvecs = MVT::Clone( *evecs, numev );
    OPT::Apply( *K, *evecs, *Kvecs );
    OPT::Apply( *M, *evecs, *Mvecs );
    MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
    // compute M-norm of residuals
    OPT::Apply( *M, *Kvecs, *Mvecs );
    MVT::MvDot( *Mvecs, *Kvecs, normV );

    os << "Number of iterations performed in LOBPCG_test.exe: " << MySolverMan.getNumIters() << std::endl
       << "Direct residual norms computed in LOBPCG_test.exe" << std::endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(M)" << std::endl
       << "----------------------------------------" << std::endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) / evals[i].realpart );
      }
      else {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) );
      }
      os << std::setw(20) << evals[i].realpart << std::setw(20) << normV[i] << std::endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (verbose && MyPID==0) {
      std::cout << std::endl << os.str() << std::endl;
    }

  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    if (verbose && MyPID==0) {
      std::cout << "End Result: TEST FAILED" << std::endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;

}
