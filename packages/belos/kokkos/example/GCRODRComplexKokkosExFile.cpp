// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a file, which must be in Matrix Market (*.mtx).  
// The problem right-hand side will be generated randomly.
//
// NOTE: No preconditioner is used in this example.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosOutputManager.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "BelosKokkosAdapter.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"
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

  typedef Kokkos::complex<double>           ST;
  typedef int                               OT;
  typedef Kokkos::DefaultExecutionSpace     EXSP;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Belos::KokkosMultiVec<ST, EXSP>         MV;
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
  int maxsubspace = 300;     // maximum number of blocks the solver can use for the subspace
  int maxrestarts = 5;       // number of restarts allowed
  bool expresidual = false;  // use explicit residual
  std::string filename("mhd1280b.mtx"); // example matrix
  std::string rhsfile("");
  MT tol = 1.0e-5;           // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("expres","impres",&expresidual,"Use explicit residual throughout.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("rhsfile",&rhsfile,"Filename for right-hand side. (*.mtx file) ");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GCRO-DR solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
  cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GCRO-DR solver.");

  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  
  // Read in a matrix Market file and use it to test the Kokkos Operator.
  KokkosSparse::CrsMatrix<ST, OT, EXSP> crsMat = 
            KokkosSparse::Impl::read_kokkos_crst_matrix<KokkosSparse::CrsMatrix<ST, OT, EXSP>>(filename.c_str());
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
  belosList.set( "Num Blocks", maxsubspace);             // Maximum number of blocks in Krylov factorization
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
  belosList.set( "Explicit Residual Test", expresidual); // Use explicit residual

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
  // ***************Start the GCRO-DR iteration*************************
  // *******************************************************************
  //
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<ST,KMV,KOP> > newSolver
    = rcp( new Belos::GCRODRSolMgr<ST,KMV,KOP,true>(rcpFromRef(problem), rcpFromRef(belosList)) );

  //
  // **********Print out information about problem*******************
  //
  std::cout << std::endl << std::endl;
  std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
  std::cout << "Number of right-hand sides: " << numrhs << std::endl;
  std::cout << "Max number of GCRO-DR iterations: " << maxiters << std::endl;
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
  std::vector<MT> actual_resids( numrhs );
  std::vector<MT> rhs_norm( numrhs );
  MV resid(numRows, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  for ( int i=0; i<numrhs; i++) {
    MT actRes = actual_resids[i]/rhs_norm[i];
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
