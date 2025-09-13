// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMaps.hpp>
#include <Galeri_XpetraMatrixTypes.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_oblackholestream.hpp>
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

// Belos
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"

// ****************************************************************************
// BEGIN RUN ROUTINE
// ****************************************************************************

template <typename ScalarType>
int run(int argc, char *argv[]) {

  // Belos solvers have the following template parameters:
  //
  //   - Scalar: The type of dot product results.
  //   - MV: The type of (multi)vectors.
  //   - OP: The type of operators (functions from multivector to
  //     multivector).  A matrix (like Tpetra::CrsMatrix) is an example
  //     of an operator; an Ifpack2 preconditioner is another example.
  //
  // Here, ST is set by the main function, MV is Tpetra::MultiVector, and OP is
  // Tpetra::Operator.
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<>::node_type;

  using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MT = typename Teuchos::ScalarTraits<ST>::magnitudeType;

  using tmap_t       = Tpetra::Map<LO,GO,NT>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  const auto comm = Tpetra::getDefaultComm();
  const int myPID = comm->getRank();

  bool verbose = false;
  bool success = true;

  auto out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  try {
    bool procVerbose = false;
    bool debug = false;
    int frequency = -1;        // frequency of status test output
    int blockSize = 1;         // blockSize
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxIters = -1;         // maximum number of iterations allowed per linear system
    int maxSubspace = 50;      // maximum number of blocks the solver can use for the subspace
    int maxRestarts = 15;      // number of restarts allowed
    int nx = 10;               // number of discretization points in each direction
    int ny = nx;
    int nz = nx;
    MT tol = 1.0e-5;           // relative residual tolerance
    std::string solverName = "Block GMRES"; // type of iterative solver

    Teuchos::CommandLineProcessor cmdp(false,true);
    Galeri::Xpetra::Parameters<GO> galeriParameters(cmdp, nx, ny, nz, "Laplace2D");
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nondebug",&debug,"Print debugging information from solver.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("block-size",&blockSize,"Block size used by GMRES.");
    cmdp.setOption("max-iters",&maxIters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    cmdp.setOption("max-subspace",&maxSubspace,"Maximum number of blocks the solver can use for the subspace.");
    cmdp.setOption("max-restarts",&maxRestarts,"Maximum number of restarts allowed for GMRES solver.");
    cmdp.setOption("solverName", &solverName, "The type of solver to use.");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    std::string watchrProblemName = std::string("Belos ") + solverName + " " + std::to_string(comm->getSize()) + " ranks";
    procVerbose = ( verbose && (myPID==0) ); // Only print on the zero processor

    if (procVerbose) {
      *out << Belos::Belos_Version() << std::endl << std::endl;
    }

    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("Driver"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

    // Set up the test problem.
    //
    // We use Trilinos' Galeri package to construct a test problem.
    // By default, we use a discretization of the 2-D Laplacian operator.
    // The global mesh size is nx * ny.
    Teuchos::ParameterList GaleriList = galeriParameters.GetParameterList();
    std::string matrixType = galeriParameters.GetMatrixType();
    RCP<tmap_t> Map;
    if (matrixType == "Laplace1D" || matrixType == "Identity") {
      Map = RCP{Galeri::Xpetra::CreateMap<LO,GO,tmap_t>("Cartesian1D", comm, GaleriList)};
    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
               matrixType == "BigStar2D" || matrixType == "AnisotropicDiffusion" || matrixType == "Elasticity2D" || matrixType == "Recirc2D") {
      Map = RCP{Galeri::Xpetra::CreateMap<LO,GO,tmap_t>("Cartesian2D", comm, GaleriList)};
    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      Map = RCP{Galeri::Xpetra::CreateMap<LO,GO,tmap_t>("Cartesian3D", comm, GaleriList)};
    }
    auto GaleriProblem = Galeri::Xpetra::BuildProblem<ST,LO,GO,tmap_t,tcrsmatrix_t,MV>(matrixType, Map, GaleriList);

    // Create matrix from problem
    auto A = GaleriProblem->BuildMatrix();

    // Create RHS using random solution vector
    RCP<MV> B = rcp (new MV (Map, numrhs));
    RCP<MV> X = rcp (new MV (Map, numrhs));
    RCP<MV> Xexact = rcp (new MV (Map, numrhs));
    MVT::MvRandom(*Xexact);

    OPT::Apply(*A, *Xexact, *B );

    // ********Other information used by block solver***********
    // *****************(can be user specified)******************

    const int numGlobalElements = B->getGlobalLength();
    if (maxIters == -1)
      maxIters = numGlobalElements/blockSize - 1; // maximum number of iterations to run

    ParameterList belosList;
    belosList.set( "Num Blocks", maxSubspace);             // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blockSize );              // BlockSize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxIters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxRestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    int verbosity = Belos::Errors + Belos::Warnings;
    if (verbose) {
      verbosity += Belos::IterationDetails + Belos::FinalSummary + Belos::StatusTestDetails;
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    if (debug) {
      verbosity += Belos::Debug;
    }
    belosList.set( "Verbosity", verbosity );

    // Construct an unpreconditioned linear problem instance.
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem();
    if (set == false) {
      if (procVerbose)
        *out << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    //
    // *******************************************************************
    // ****************Start the solver iteration*************************
    // *******************************************************************
    //
    // Create a solver factory
    Belos::SolverFactory<double,MV,OP> factory;

    // Create an iterative solver manager
    RCP< Belos::SolverManager<double,MV,OP> > newSolver = factory.create (solverName, rcp(&belosList,false));

    // Set the problem on the solver manager
    newSolver->setProblem( rcp(&problem,false) );

    // **********Print out information about problem*******************

    if (procVerbose) {
      *out << std::endl << std::endl;
      *out << "Solver: " << solverName << std::endl;
      *out << "Dimension of matrix: " << numGlobalElements << std::endl;
      *out << "Number of right-hand sides: " << numrhs << std::endl;
      *out << "Block size used by solver: " << blockSize << std::endl;
      *out << "Max number of restarts allowed: " << maxRestarts << std::endl;
      *out << "Max number of iterations per linear system: " << maxIters << std::endl;
      *out << "Relative residual tolerance: " << tol << std::endl;
      *out << std::endl;
    }

    // Perform solve
    Belos::ReturnType ret = newSolver->solve();

    stacked_timer->stopBaseTimer();
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(*out, comm, options);
    auto xmlOut = stacked_timer->reportWatchrXML(watchrProblemName, comm);
    if (xmlOut.length())
      *out << "\nAlso created Watchr performance report " << xmlOut << '\n';

    // Get the number of iterations for this solve.
    int numIters = newSolver->getNumIters();
    if (procVerbose)
      *out << "Number of iterations performed for this solve: " << numIters << std::endl;

    // Compute actual residuals.
    bool badRes = false;
    std::vector<ST> actualResids( numrhs );
    std::vector<ST> rhsNorm( numrhs );
    MV resid(Map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actualResids );
    MVT::MvNorm( *B, rhsNorm );
    if (procVerbose) {
      *out<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        *out<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret!=Belos::Converged || badRes) {
      success = false;
      if (procVerbose)
        *out << "End Result: TEST FAILED" << std::endl;
    } else {
      if (procVerbose)
        *out << "End Result: TEST PASSED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[]) {
  // run with different ST
  return run<double>(argc,argv);
  // run<float>(argc,argv); // FAILS
}
