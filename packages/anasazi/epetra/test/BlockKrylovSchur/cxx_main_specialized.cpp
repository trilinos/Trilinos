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
// This test is for BlockKrylovSchur solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the BlockKrylovSchurSolMgr solver manager.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziSpecializedEpetraAdapter.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBasicSort.hpp"
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
  bool shortrun = false;
  bool insitu = false;
  std::string which("LM");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (debug) verbose = true;

  typedef double ScalarType;
  typedef ScalarTraits<ScalarType>                     SCT;
  typedef SCT::magnitudeType                           MagnitudeType;
  typedef Anasazi::MultiVec<ScalarType>                MV;
  typedef Anasazi::Operator<ScalarType>                OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>       MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>    OPT;

  const ScalarType ONE  = SCT::one();

  if (verbose && MyPID == 0) {
    std::cout << Anasazi::Anasazi_Version() << std::endl << std::endl;
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
  const int maxIterCG = 100;
  const double tolCG = 1e-7;
  
  RCP<BlockPCGSolver> opStiffness = rcp( new BlockPCGSolver(Comm, M.get(), tolCG, maxIterCG, 0) );
  opStiffness->setPreconditioner( 0 );
  RCP<const Anasazi::EpetraGenOp> InverseOp = rcp( new Anasazi::EpetraGenOp( opStiffness, K ) );

  // Create the initial vectors
  int blockSize = 3;
  RCP<Anasazi::EpetraOpMultiVec> ivec = rcp( new Anasazi::EpetraOpMultiVec(K, K->OperatorDomainMap(), blockSize) );
  ivec->MvRandom();

  // Create eigenproblem
  const int nev = 5;
  RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(InverseOp, ivec) );
  //
  // Inform the eigenproblem that the operator InverseOp is Hermitian under an M inner-product
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
  int numBlocks;
  int maxRestarts;
  if (shortrun) {
    maxRestarts = 25;
    numBlocks = 5;
  }
  else {
    maxRestarts = 50;
    numBlocks = 10;
  }
  int stepSize = numBlocks*maxRestarts;
  MagnitudeType tol = tolCG * 10.0;
  // Create a sort manager to pass into the block Krylov-Schur solver manager
  // -->  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
  // -->  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
  //      so you can also pass in the parameter "Which", instead of a sort manager.
  RCP<Anasazi::SortManager<MagnitudeType> > MySort =     
    rcp( new Anasazi::BasicSort<MagnitudeType>( which ) );
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Sort Manager", MySort );
  //MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Step Size", stepSize );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "In Situ Restarting", insitu );
  //
  // Create the solver manager
  Anasazi::BlockKrylovSchurSolMgr<ScalarType,MV,OP> MySolverMgr(problem, MyPL);
  // 
  // Check that the parameters were all consumed
  if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
      MyPL.getEntryPtr("Block Size")->isUsed() == false ||
      MyPL.getEntryPtr("Num Blocks")->isUsed() == false ||
      MyPL.getEntryPtr("Maximum Restarts")->isUsed() == false ||
      MyPL.getEntryPtr("Step Size")->isUsed() == false ||
      MyPL.getEntryPtr("In Situ Restarting")->isUsed() == false ||
      MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false) {
    if (verbose && MyPID==0) {
      std::cout << "Failure! Unused parameters: " << std::endl;
      MyPL.unused(std::cout);
    }
  }


  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged && shortrun==false) {
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

    RCP<OP> KOp = rcp( new Anasazi::EpetraOp( K ));    
    RCP<OP> MOp = rcp( new Anasazi::EpetraOp( M ));    
    RCP<MV> Mvecs = MVT::Clone( *evecs, numev ),
                    Kvecs = MVT::Clone( *evecs, numev );
    OPT::Apply( *KOp, *evecs, *Kvecs );
    OPT::Apply( *MOp, *evecs, *Mvecs );
    MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
    // compute 2-norm of residuals
    std::vector<MagnitudeType> resnorm(numev);
    MVT::MvNorm( *Kvecs, resnorm );

    os << "Number of iterations performed in BlockKrylovSchur_test.exe: " << MySolverMgr.getNumIters() << std::endl
       << "Direct residual norms computed in BlockKrylovSchur_test.exe" << std::endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual" << std::endl
       << "----------------------------------------" << std::endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(evals[i].realpart) != SCT::zero() ) {
        resnorm[i] = SCT::magnitude( SCT::squareroot( resnorm[i] ) / evals[i].realpart );
      }
      else {
        resnorm[i] = SCT::magnitude( SCT::squareroot( resnorm[i] ) );
      }
      os << std::setw(20) << evals[i].realpart << std::setw(20) << resnorm[i] << std::endl;
      if ( resnorm[i] > tol ) {
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
