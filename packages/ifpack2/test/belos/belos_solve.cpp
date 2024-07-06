// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_Comm.hpp"

#include "Ifpack2_Parameters.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers_decl.hpp"

#include "build_problem.hpp"
#include "build_solver.hpp"

void
process_command_line (bool& printedHelp,
                      std::string& xml_file,
                      bool& useStackedTimer,
                      std::string& problem_name,
                      int argc,
                      char*argv[])
{
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("xml_file", &xml_file, "XML Parameters file");
  cmdp.setOption("with_stacked_timer", "without_stacked_timer", &useStackedTimer,
      "Whether to run with a StackedTimer and print the timer tree at the end");
  cmdp.setOption("problem_name", &problem_name, "Problem name for Watchr plot");

  const auto result = cmdp.parse (argc, argv);

  // mfh 21 Apr 2016: By ignoring options that this executable doesn't
  // recognize, we can pass them through to (e.g.,) Kokkos.

  if (result == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    printedHelp = true; // not an error to ask for help
  }
  else if (result == Teuchos::CommandLineProcessor::PARSE_ERROR) {
    throw std::runtime_error ("Error parsing command-line.");
  }
}

int main (int argc, char* argv[])
{
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>        STS;
  typedef STS::magnitudeType                   Magnitude;
  typedef Tpetra::Map<>::local_ordinal_type    LO;
  typedef Tpetra::Map<>::global_ordinal_type   GO;
  typedef Tpetra::Map<>::node_type             Node;
  typedef Tpetra::MultiVector<Scalar,LO,GO>    TMV;
  typedef Tpetra::Operator<Scalar,LO,GO>       TOP;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
  typedef Belos::LinearProblem<Scalar,TMV,TOP> BLinProb;
  typedef Belos::SolverManager<Scalar,TMV,TOP> BSolverMgr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::StackedTimer;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  bool success = true;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    Teuchos::Time timer("total");
    timer.start();

    RCP<const Teuchos::Comm<int> > comm =
      Tpetra::getDefaultComm();


    //Just get one parameter from the command-line: the name of an xml file
    //to get parameters from.

    std::string xml_file("calore1_mm.xml");
    std::string problem_name("Belos-Ifpack2 Setup-Solve");
    bool useStackedTimer = false;
    {
      bool printedHelp = false;
      process_command_line (printedHelp, xml_file, useStackedTimer, problem_name, argc, argv);
      if (printedHelp) {
        return EXIT_SUCCESS;
      }
    }

    //Read the contents of the xml file into a ParameterList. That parameter list
    //should specify a matrix-file and optionally which Belos solver to use, and
    //which Ifpack2 preconditioner to use, etc. If there are sublists of parameters
    //for Belos and Ifpack2, those will be passed to the respective destinations
    //from within the build_problem and build_solver functions.

    *out << "Every proc reading parameters from xml_file: "
         << xml_file << std::endl;
    Teuchos::ParameterList test_params =
      Teuchos::ParameterXMLFileReader(xml_file).getParameters();

    //The build_problem function is located in build_problem.hpp.
    //Note that build_problem calls build_precond and sets a preconditioner on the
    //linear-problem, if a preconditioner is specified.

    RCP<BLinProb> problem =
      build_problem<Scalar,LO,GO,Node>(test_params, comm);

    RCP<StackedTimer> stackedTimer;
    if(useStackedTimer)
    {
      stackedTimer = rcp(new StackedTimer("Belos-Ifpack2 Setup-Solve"));
      Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
    }

    // defaullt convergence tol
    Magnitude tol = 1e-7;
    Teuchos::ParameterList& BelosParams = test_params.sublist("Belos");
    if (BelosParams.isParameter ("Convergence Tolerance")) {
      tol = BelosParams.get<Magnitude> ("Convergence Tolerance");
    } else {
      BelosParams.set<Magnitude> ("Convergence Tolerance", tol);
    }

    //Build the preconditioner (if one is enabled), and bind it to the problem
    std::string tifpack_precond("not specified");
    Ifpack2::getParameter (test_params, "Ifpack2::Preconditioner", tifpack_precond);
    std::string prec_side("Left");
    Ifpack2::getParameter (test_params, "Preconditioner Side", prec_side);
    if (tifpack_precond != "not specified") {
      RCP<TOP> precond;
      if (tifpack_precond == "RBILUK") {
        int blockSize = 0;
        Teuchos::ParameterList& prec_params = test_params.sublist("Ifpack2");
        Ifpack2::getParameter (prec_params, "fact: block size", blockSize);
        assert(blockSize >= 1);
        auto A_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(problem->getOperator());
        if (blockSize > 1) {
          auto crs_matrix_block_filled = Tpetra::fillLogicalBlocks(*A_crs, blockSize);
          auto A = Teuchos::rcp_const_cast<const block_crs_matrix_type>(Tpetra::convertToBlockCrsMatrix(*crs_matrix_block_filled, blockSize));
          precond = build_precond<Scalar,LO,GO,Node> (test_params, A);
        }
        else {
          auto A = Teuchos::rcp_const_cast<const block_crs_matrix_type>(Tpetra::convertToBlockCrsMatrix(*A_crs, blockSize));
          precond = build_precond<Scalar,LO,GO,Node> (test_params, A);
        }
      }
      else {
        auto A = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(problem->getOperator());
        precond = build_precond<Scalar,LO,GO,Node> (test_params, A);
      }
      if (prec_side == "Left")
        problem->setLeftPrec (precond);
      else if (prec_side == "Right")
        problem->setRightPrec (precond);
    }

    problem->setProblem ();

    //The build_solver function is located in build_solver.hpp:

    RCP<BSolverMgr> solver = build_solver<Scalar,TMV,TOP>(test_params, problem);

    Belos::ReturnType ret = solver->solve();

    *out << "Converged in " << solver->getNumIters() << " iterations." << std::endl;

    RCP<const TOP> prec = problem->getLeftPrec();
    if (prec !=Teuchos::null) {
      *out << "Preconditioner attributes:" << std::endl;
      prec->describe (*out, Teuchos::VERB_LOW);
    }
    prec = problem->getRightPrec();
    if (prec !=Teuchos::null) {
      *out << "Preconditioner attributes:" << std::endl;
      prec->describe (*out, Teuchos::VERB_LOW);
    }

    RCP<TMV> R = rcp(new TMV(*problem->getRHS(), Teuchos::Copy));
    problem->computeCurrResVec(&*R, &*problem->getLHS(), &*problem->getRHS());
    Teuchos::Array<Magnitude> normsR(R->getNumVectors());
    Teuchos::Array<Magnitude> normsB(R->getNumVectors());
    R->norm2(normsR);
    problem->getRHS()->norm2(normsB);

    if (normsR.size() < 1) {
      throw std::runtime_error("ERROR: norms.size()==0 indicates R->getNumVectors()==0.");
    }

    *out << "2-Norm of 0th RHS      vec: " << normsB[0] << std::endl;
    *out << "2-Norm of 0th residual vec: " << normsR[0] << " -> " << normsR[0]/normsB[0] << std::endl;
    *out << "Achieved tolerance: " << solver->achievedTol() << ", Requested tolerance: " << tol << std::endl;
    normsR[0] /= normsB[0];

    //If the xml file specified a number of iterations to expect, then we will
    //use that as a test pass/fail criteria.

    if (test_params.isParameter("expectNumIters")) {
      int expected_iters = 0;
      Ifpack2::getParameter(test_params, "expectNumIters", expected_iters);
      int actual_iters = solver->getNumIters();
      if (ret == Belos::Converged && actual_iters <= expected_iters && normsR[0] < tol) {
      }
      else {
        success = false;
        *out << "Actual iters(" << actual_iters << ") > expected number of iterations ("
             << expected_iters <<"), or resid-norm(" << normsR[0] << ") >= " << tol << std::endl;
      }
    }

    timer.stop();
    *out << "proc 0 total program time: " << timer.totalElapsedTime()
         << std::endl;

    if(useStackedTimer)
    {
      stackedTimer->stopBaseTimer();
      StackedTimer::OutputOptions options;
      options.num_histogram=3;
      options.print_warnings = false;
      options.output_histogram = true;
      options.output_fraction=true;
      options.output_minmax = true;
      stackedTimer->report(std::cout, comm, options);
      auto xmlOut = stackedTimer->reportWatchrXML(problem_name, comm);
      if(comm->getRank() == 0)
      {
        if(xmlOut.length())
          std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success)

  if (success) {
    *out << "End Result: TEST PASSED\n";
  }
  else {
    *out << "End Result: TEST FAILED\n";
  }

  return ( success ? 0 : 1 );
}

