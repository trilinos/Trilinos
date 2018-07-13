// @HEADER
// ***********************************************************************
//
//                 Belos: Block Eigensolvers Package
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

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

// The Trilinos package Galeri has many example problems.
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Include selected communicator class required by Epetra objects
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // EPETRA_MPI

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int
main (int argc, char *argv[])
{
  int MyPID = 0;

  // Belos solvers have the following template parameters:
  //
  //   - Scalar: The type of dot product results.
  //   - MV: The type of (multi)vectors.
  //   - OP: The type of operators (functions from multivector to
  //     multivector).  A matrix (like Epetra_CrsMatrix) is an example
  //     of an operator; an Ifpack preconditioner is another example.
  //
  // Here, Scalar is double, MV is Epetra_MultiVector, and OP is
  // Epetra_Operator.
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef EPETRA_MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif // EPETRA_MPI

bool verbose = false;
bool success = true;
try {
  bool proc_verbose = false;
  bool debug = false;
  bool pseudo = false;       // use pseudo-block or block gmres
  int frequency = -1;        // frequency of status test output
  int blocksize = 1;         // blocksize
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  int maxsubspace = 50;      // maximum number of blocks the solver can use for the subspace
  int maxrestarts = 15;      // number of restarts allowed
  int nx = 10;               // number of discretization points in each direction
  MT tol = 1.0e-5;           // relative residual tolerance
  std::string ortho = "DGKS";// orthogonalization method

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nondebug",&debug,"Print debugging information from solver.");
  cmdp.setOption("pseudo","block",&pseudo,"Use pseudo-block or block GMRES solver.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size used by GMRES.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
  cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
  cmdp.setOption("nx",&nx,"Number of discretization points in each direction of 3D Laplacian.");
  cmdp.setOption("ortho",&ortho,"Orthogonalization being used by GMRES solver.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose

  //
  // Set up the test problem.
  //
  // We use Trilinos' Galeri package to construct a test problem.
  // Here, we use a discretization of the 2-D Laplacian operator.
  // The global mesh size is nx * nx.
  //
  Teuchos::ParameterList GaleriList;
  GaleriList.set ("n", nx * nx * nx);
  GaleriList.set ("nx", nx);
  GaleriList.set ("ny", nx);
  GaleriList.set ("nz", nx);
  RCP<Epetra_Map> Map = rcp (Galeri::CreateMap ("Linear", Comm, GaleriList));
  RCP<Epetra_RowMatrix> A =
    rcp (Galeri::CreateCrsMatrix ("Laplace3D", &*Map, GaleriList));

  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

  // Create RHS using random solution vector
  RCP<MV> B = rcp (new MV (*Map, numrhs));
  RCP<MV> X = rcp (new MV (*Map, numrhs));
  RCP<MV> Xexact = rcp (new MV (*Map, numrhs));
  Xexact->Random ();

  A->Apply( *Xexact, *B );

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Num Blocks", maxsubspace);             // Maximum number of blocks in Krylov factorization
  belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Orthogonalization", ortho );           // Orthogonalization used by iterative solver
  int verbosity = Belos::Errors + Belos::Warnings;
  if (verbose) {
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  if (debug) {
    verbosity += Belos::Debug;
  }
  belosList.set( "Verbosity", verbosity );
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<double,MV,OP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  //
  //
  // *******************************************************************
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<double,MV,OP> > newSolver;
  if (pseudo && (blocksize == 1))
    newSolver = rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));
  else
    newSolver = rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of restarts allowed: " << maxrestarts << std::endl;
    std::cout << "Max number of Gmres iterations per linear system: " << maxiters << std::endl;
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = newSolver->solve();
  //
  // Get the number of iterations for this solve.
  //
  int numIters = newSolver->getNumIters();
  if (proc_verbose)
    std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;
  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(*Map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      if (actRes > tol) badRes = true;
    }
  }

if (ret!=Belos::Converged || badRes) {
  success = false;
  if (proc_verbose)
    std::cout << "End Result: TEST FAILED" << std::endl;
} else {
  if (proc_verbose)
    std::cout << "End Result: TEST PASSED" << std::endl;
}
}
TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
MPI_Finalize();
#endif

return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
