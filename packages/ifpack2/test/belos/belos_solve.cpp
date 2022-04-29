/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//@HEADER
*/

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
  typedef Tpetra::Map<>::local_ordinal_type    LO;
  typedef Tpetra::Map<>::global_ordinal_type   GO;
  typedef Tpetra::Map<>::node_type             Node;
  typedef Tpetra::MultiVector<Scalar,LO,GO>    TMV;
  typedef Tpetra::Operator<Scalar,LO,GO>       TOP;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
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

    //Build the preconditioner (if one is enabled), and bind it to the problem
    std::string tifpack_precond("not specified");
    Ifpack2::getParameter (test_params, "Ifpack2::Preconditioner", tifpack_precond);
    std::string prec_side("Left");
    Ifpack2::getParameter (test_params, "Preconditioner Side", prec_side);
    if (tifpack_precond != "not specified") {
      auto A = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(problem->getOperator());
      RCP<TOP> precond = build_precond<Scalar,LO,GO,Node> (test_params, A);
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

    RCP<TMV> R = rcp(new TMV(*problem->getRHS()));
    problem->computeCurrResVec(&*R, &*problem->getLHS(), &*problem->getRHS());
    Teuchos::Array<Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(R->getNumVectors());
    R->norm2(norms);

    if (norms.size() < 1) {
      throw std::runtime_error("ERROR: norms.size()==0 indicates R->getNumVectors()==0.");
    }

    *out << "2-Norm of 0th residual vec: " << norms[0] << std::endl;
    *out << "Achieved tolerance: " << solver->achievedTol() << std::endl;

    //If the xml file specified a number of iterations to expect, then we will
    //use that as a test pass/fail criteria.

    if (test_params.isParameter("expectNumIters")) {
      int expected_iters = 0;
      Ifpack2::getParameter(test_params, "expectNumIters", expected_iters);
      int actual_iters = solver->getNumIters();
      if (ret == Belos::Converged && actual_iters <= expected_iters && norms[0] < 1.e-7) {
      }
      else {
        success = false;
        *out << "Actual iters("<<actual_iters
             <<") > expected number of iterations ("
             <<expected_iters<<"), or resid-norm(" << norms[0] << ") >= 1.e-7"<<std::endl;
      }
    }

    timer.stop();
    *out << "proc 0 total program time: " << timer.totalElapsedTime()
         << std::endl;

    if(useStackedTimer)
    {
      stackedTimer->stopBaseTimer();
      StackedTimer::OutputOptions options;
      options.print_warnings = false;
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

