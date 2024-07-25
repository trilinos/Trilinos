// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Test for BlockGCRODRSolMgr (Block Recycling GMRES, by Kirk
// Soodhalter and Michael Parks).  This just tests compilation and
// setParameters() for now.  Later, we'll test actually solving linear
// systems.
//

#include "BelosConfigDefs.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosTpetraTestFramework.hpp"
#include "BelosBlockGCRODRSolMgr.hpp"

#include <Tpetra_Core.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

template <typename ScalarType>
int run (int argc, char *argv[])
{
  // Teuchos usings
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::oblackholestream;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  
  // Typedefs for Tpetra template arguments.
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<ST, LO, GO>::node_type;
  
  // Tpetra objects which are the MV and OP template parameters of the
  // Belos specialization which we are testing.
  using MV = typename Tpetra::MultiVector<ST, LO, GO>;
  using OP = typename Tpetra::Operator<ST, LO, GO>;
  
  // Other typedefs.
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST, LO, GO>;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<oblackholestream> blackHole (new oblackholestream);
  const int myRank = comm->getRank ();

  // Output stream that prints only on Process 0.
  RCP<FancyOStream> out;
  if (myRank == 0) {
    out = Teuchos::getFancyOStream (rcpFromRef (std::cout));
  } else {
    out = Teuchos::getFancyOStream (blackHole);
  }

  // Get test parameters from command-line processor.
  // CommandLineProcessor always understands int, but may not
  // understand GO.  We convert to the latter below.
  int numRows = comm->getSize() * 100;
  bool tolerant = false;
  bool verbose = false;
  bool debug = false;
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption("numRows", &numRows,
                 "Global number of rows (and columns) in the sparse matrix to generate.");
  cmdp.setOption("tolerant", "intolerant", &tolerant,
                 "Whether to parse files tolerantly.");
  cmdp.setOption("verbose", "quiet", &verbose,
                 "Print messages and results.");
  cmdp.setOption("debug", "release", &debug,
                 "Run debugging checks and print copious debugging output.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    *out << "\nEnd Result: TEST FAILED" << std::endl;
    return EXIT_FAILURE;
  }
  // Output stream for verbose output.
  RCP<FancyOStream> verbOut = verbose ? out : getFancyOStream (blackHole);

  const bool success = true;

  // Test whether it's possible to instantiate the solver.
  // This is a minimal compilation test.
  *verbOut << "Instantiating Block GCRODR solver" << std::endl;
  Belos::BlockGCRODRSolMgr<ST, MV, OP> solver;

  // Test setting solver parameters. For now, we just use an empty
  // (but non-null) parameter list, which the solver should fill in
  // with defaults.
  *verbOut << "Setting solver parameters" << std::endl;
  RCP<ParameterList> solverParams = parameterList ();
  solver.setParameters (solverParams);
  
  // Create a linear system to solve.
  *verbOut << "Creating linear system" << std::endl;
  RCP<tcrsmatrix_t> A;
  RCP<MV> X_guess, X_exact, B;
  {
    Teuchos::RCP<NT> node; // can be null; only for type deduction
    Belos::Tpetra::ProblemMaker<tcrsmatrix_t> factory (comm, node, out, tolerant, debug);

    RCP<ParameterList> problemParams = parameterList ();
    problemParams->set ("Global number of rows",
                        static_cast<GO> (numRows));
    problemParams->set ("Problem type", std::string ("Nonsymmetric"));
    factory.makeProblem (A, X_guess, X_exact, B, problemParams);
  }

  // Approximate solution vector is a copy of the guess vector.
  RCP<MV> X (new MV (*X_guess));

  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::logic_error,
                             "The sparse matrix is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(X_guess.is_null(), std::logic_error,
                             "The initial guess X_guess is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(X_exact.is_null(), std::logic_error,
                             "The exact solution X_exact is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.is_null(), std::logic_error,
                             "The right-hand side B is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(X.is_null(), std::logic_error,
                             "The approximate solution vector X is null!");

  typedef Belos::LinearProblem<ST, MV, OP> problem_type;
  RCP<problem_type> problem (new problem_type (A, X, B));
  problem->setProblem ();
  solver.setProblem (problem);

  *verbOut << "Solving linear system" << std::endl;
  Belos::ReturnType result = solver.solve ();

  *verbOut << "Result of solve: "
           << Belos::convertReturnTypeToString (result)
           << std::endl;

  // Make sure that all the processes finished.
  comm->barrier ();

  if (success) {
    *out << "\nEnd Result: TEST PASSED" << std::endl;
    return EXIT_SUCCESS;
  }
  else {
    *out << "\nEnd Result: TEST FAILED" << std::endl;
    return EXIT_FAILURE;
  }
}

int main(int argc, char *argv[]) {
  return run<double>(argc, argv);
  // return run<float>(argc, argv);
}

