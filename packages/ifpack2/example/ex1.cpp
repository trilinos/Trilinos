// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <Ifpack2_Factory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_TimeMonitor.hpp>

namespace { // (anonymous)

// Values of command-line arguments.
struct CmdLineArgs {
  CmdLineArgs () :
    usePreconditioner (true)
  {}

  std::string matrixFilename;
  std::string rhsFilename;
  bool usePreconditioner;
};

// Read in values of command-line arguments.
bool
getCmdLineArgs (CmdLineArgs& args, int argc, char* argv[])
{
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("usePreconditioner", "noPreconditioner",
                  &args.usePreconditioner, "Whether to use a preconditioner");
  cmdp.setOption ("matrixFilename", &args.matrixFilename, "Name of Matrix "
                  "Market file with the sparse matrix A");
  cmdp.setOption ("rhsFilename", &args.rhsFilename, "Name of Matrix Market "
                  "file with the right-hand side vector(s) B");
  auto result = cmdp.parse (argc, argv);
  return result == Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL;
}

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using std::cerr;
  using std::endl;
  typedef Tpetra::CrsMatrix<> crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::MultiVector<> MV;
  typedef Tpetra::Operator<> OP;
  typedef Tpetra::RowMatrix<> row_matrix_type;
  typedef MV::scalar_type scalar_type;
  typedef Ifpack2::Preconditioner<> prec_type;
  typedef Belos::LinearProblem<scalar_type, MV, OP> problem_type;
  typedef Belos::SolverManager<scalar_type, MV, OP> solver_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  RCP<Time> totalTime = Teuchos::TimeMonitor::getNewTimer ("Total");
  RCP<Time> precSetupTime =
    Teuchos::TimeMonitor::getNewTimer ("Preconditioner setup");
  RCP<Time> probSetupTime =
    Teuchos::TimeMonitor::getNewTimer ("Problem setup");
  RCP<Time> solveTime = Teuchos::TimeMonitor::getNewTimer ("Solve");

  Teuchos::TimeMonitor totalTimeMon (*totalTime);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  // Get command-line arguments.
  CmdLineArgs args;
  const bool gotCmdLineArgs = getCmdLineArgs (args, argc, argv);
  if (! gotCmdLineArgs) {
    if (comm->getRank () == 0) {
      cerr << "Failed to get command-line arguments!" << endl;
    }
    return EXIT_FAILURE;
  }
  if (args.matrixFilename == "") {
    if (comm->getRank () == 0) {
      cerr << "Must specify sparse matrix filename!" << endl;
    }
    return EXIT_FAILURE;
  }
  if (args.rhsFilename == "") {
    if (comm->getRank () == 0) {
      cerr << "Must specify filename for loading right-hand side(s)!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Read sparse matrix A from Matrix Market file.
  RCP<crs_matrix_type> A =
    reader_type::readSparseFile (args.matrixFilename, comm);
  if (A.is_null ()) {
    if (comm->getRank () == 0) {
      cerr << "Failed to load sparse matrix A from file "
        "\"" << args.matrixFilename << "\"!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Read right-hand side vector(s) B from Matrix Market file.
  RCP<const map_type> map = A->getRangeMap ();
  RCP<MV> B = reader_type::readDenseFile (args.rhsFilename, comm, map);
  if (B.is_null ()) {
    if (comm->getRank () == 0) {
      cerr << "Failed to load right-hand side vector(s) from file \""
           << args.rhsFilename << "\"!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Create Belos iterative linear solver.
  RCP<solver_type> solver;
  RCP<ParameterList> solverParams (new ParameterList ());
  {
    Belos::SolverFactory<scalar_type, MV, OP> belosFactory;
    solver = belosFactory.create ("GMRES", solverParams);
  }
  if (solver.is_null ()) {
    if (comm->getRank () == 0) {
      cerr << "Failed to create Belos solver!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Optionally, create Ifpack2 preconditioner.
  RCP<prec_type> M;
  if (args.usePreconditioner) {
    Teuchos::TimeMonitor precSetupTimeMon (*precSetupTime);
    M = Ifpack2::Factory::create<row_matrix_type> ("RELAXATION", A);
    if (M.is_null ()) {
      if (comm->getRank () == 0) {
        cerr << "Failed to create Ifpack2 preconditioner!" << endl;
      }
      return EXIT_FAILURE;
    }
    M->initialize ();
    M->compute ();
  }

  // Set up the linear problem to solve.
  RCP<MV> X (new MV (A->getDomainMap (), B->getNumVectors ()));
  RCP<problem_type> problem;
  {
    Teuchos::TimeMonitor probSetupTimeMon (*probSetupTime);
    problem = rcp (new problem_type (A, X, B));
    if (! M.is_null ()) {
      problem->setRightPrec (M);
    }
    problem->setProblem ();
    solver->setProblem (problem);
  }

  // Solve the linear system.
  {
    Teuchos::TimeMonitor solveTimeMon (*solveTime);
    Belos::ReturnType solveResult = solver->solve ();
    if (solveResult != Belos::Converged) {
      if (comm->getRank () == 0) {
        cerr << "Solve failed to converge!" << endl;
      }
      return EXIT_FAILURE;
    }
  }

  // Report timings.
  Teuchos::TimeMonitor::report (comm.ptr (), std::cout);

  return EXIT_SUCCESS;
}
