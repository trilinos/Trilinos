//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
//@HEADER
//
// This driver reads a problem from a file, which must be in Matrix Market (*.mtx).  
// The problem right-hand side will be generated randomly.
//
// NOTE: No preconditioner is used in this example.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosOutputManager.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "BelosKokkosAdapter.hpp"
#include "KokkosKernels_IOUtils.hpp"
#ifdef HAVE_MPI
  #include <mpi.h>
#endif

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

bool success = true;
  Kokkos::initialize();
  {

  typedef double                            ST;
  typedef int                               OT;
  typedef Kokkos::DefaultExecutionSpace     EXSP;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Belos::KokkosMultiVec<ST, EXSP>         MV;
  typedef Belos::KokkosCrsOperator<ST, OT, EXSP>       OP;
  typedef Belos::MultiVec<ST> KMV;
  typedef Belos::Operator<ST> KOP; 
  typedef Belos::MultiVecTraits<ST,KMV>     MVT;
  typedef Belos::OperatorTraits<ST,KMV,KOP>  OPT;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

bool verbose = true;
try {
  int frequency = 25;        // frequency of status test output.
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  bool expresidual = false; // use explicit residual
  std::string filename("bcsstk12.mtx"); // example matrix
  MT tol = 1.0e-5;           // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("expres","impres",&expresidual,"Use explicit residual throughout.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by Gmres solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");

  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  
  // Read in a matrix Market file and use it to test the Kokkos Operator.
  KokkosSparse::CrsMatrix<ST, OT, EXSP> crsMat = 
            KokkosKernels::Impl::read_kokkos_crst_matrix<KokkosSparse::CrsMatrix<ST, OT, EXSP>>(filename.c_str()); 
  RCP<Belos::KokkosCrsOperator<ST, OT, EXSP>> A = 
            rcp(new Belos::KokkosCrsOperator<ST,OT,EXSP>(crsMat));
  OT numRows = crsMat.numRows();

  Teuchos::RCP<MV> X = Teuchos::rcp( new MV(numRows, numrhs) );
  X->MvRandom();
  Teuchos::RCP<MV> B = Teuchos::rcp( new MV(numRows, numrhs) );
  OPT::Apply(*A,*X,*B);
  X->MvInit(0.0);

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GetGlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements - 1; // maximum number of iterations to run
  
  ParameterList belosList;
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Explicit Residual Test", expresidual);      // use explicit residual

  if (verbose) {
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
		   Belos::StatusTestDetails + Belos::FinalSummary + Belos::TimingDetails);
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  else
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
  
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<ST,KMV,KOP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  //
  // *******************************************************************
  // **************Start the block CG iteration*************************
  // *******************************************************************
  //
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<ST,KMV,KOP> > newSolver
    = rcp( new Belos::BlockCGSolMgr<ST,KMV,KOP>(rcpFromRef(problem), rcpFromRef(belosList)) );

  //
  // **********Print out information about problem*******************
  //
  std::cout << std::endl << std::endl;
  std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
  std::cout << "Number of right-hand sides: " << numrhs << std::endl;
  std::cout << "Max number of Gmres iterations: " << maxiters << std::endl;
  std::cout << "Relative residual tolerance: " << tol << std::endl;
  std::cout << std::endl;
  //
  // Perform solve
  //
  Belos::ReturnType ret;
  ret = newSolver->solve();

  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<ST> actual_resids( numrhs );
  std::vector<ST> rhs_norm( numrhs );
  MV resid(numRows, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  for ( int i=0; i<numrhs; i++) {
    ST actRes = actual_resids[i]/rhs_norm[i];
    std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
    if (actRes > tol) badRes = true;
  }

  if (ret!=Belos::Converged || badRes) {
    success = false;
    std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
  } else {
    success = true;
    std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
  }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
