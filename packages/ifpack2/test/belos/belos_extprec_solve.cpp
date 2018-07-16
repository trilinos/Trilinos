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

#include "Ifpack2_ConfigDefs.hpp"

#ifdef HAVE_IFPACK2_QD

#include "Tpetra_Core.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#include "Ifpack2_Parameters.hpp"

#include <qd/dd_real.h>
// Set a flag that we are using QD so that the solver builder only instantiates supported solvers
#define USING_QD

#include "build_problem.hpp"
#include "build_solver.hpp"

void process_command_line(int argc, char*argv[], std::string& xml_file);

int main(int argc, char*argv[])
{
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  Teuchos::Time timer("total");
  timer.start();

  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  typedef dd_real Scalar;
  typedef int LO; //LocalOrdinal
  typedef int GO; //GlobalOrdinal
  typedef Tpetra::Map<LO, GO>::node_type Node;
  typedef Tpetra::MultiVector<Scalar,LO,GO> TMV;
  typedef Tpetra::Operator<Scalar,LO,GO>    TOP;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>   BLinProb;
  typedef Belos::SolverManager<Scalar,TMV,TOP>   BSolverMgr;

  //Just get one parameter from the command-line: the name of an xml file
  //to get parameters from.

  std::string xml_file("calore1_mm.xml");
  process_command_line(argc, argv, xml_file);

  //Read the contents of the xml file into a ParameterList. That parameter list
  //should specify a matrix-file and optionally which Belos solver to use, and
  //which Ifpack2 preconditioner to use, etc. If there are sublists of parameters
  //for Belos and Ifpack2, those will be passed to the respective destinations
  //from within the build_problem and build_solver functions.

  std::cout << "Every proc reading parameters from xml_file: "
            << xml_file << std::endl;
  Teuchos::ParameterList test_params =
      Teuchos::ParameterXMLFileReader(xml_file).getParameters();

  //The build_problem function is located in build_problem.hpp.
  //Note that build_problem calls build_precond and sets a preconditioner on the
  //linear-problem, if a preconditioner is specified.

  Teuchos::RCP<BLinProb> problem = build_problem<Scalar,LO,GO,Node>(test_params, comm);

  //The build_solver function is located in build_solver.hpp:

  Teuchos::RCP<BSolverMgr> solver = build_solver<Scalar,TMV,TOP>(test_params, problem);

  Belos::ReturnType ret = solver->solve();

  if (comm->getRank() == 0) {
    std::cout << "Converged in " << solver->getNumIters() << " iterations." << std::endl;
  }

  //Next compute residual vector and then 2-norm of residual:

  Teuchos::RCP<TMV> R = Teuchos::rcp(new TMV(*problem->getRHS()));
  problem->computeCurrResVec(&*R, &*problem->getLHS(), &*problem->getRHS());
  Teuchos::Array<Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(R->getNumVectors());
  R->norm2(norms);

  if (norms.size() < 1) {
    throw std::runtime_error("ERROR: norms.size()==0 indicates R->getNumVectors()==0.");
  }

  if (comm->getRank() == 0) {
    std::cout << "2-Norm of 0th residual vec: " << norms[0] << std::endl;
  }

  //If the xml file specified a number of iterations to expect, then we will
  //use that as a test pass/fail criteria.

  if (test_params.isParameter("expectNumIters")) {
    int expected_iters = 0;
    Ifpack2::getParameter(test_params, "expectNumIters", expected_iters);
    int actual_iters = solver->getNumIters();
    if (ret == Belos::Converged && actual_iters <= expected_iters && norms[0] < 1.e-7) {
      if (comm->getRank() == 0) {
        std::cout << "End Result: TEST PASSED" << std::endl;
      }
    }
    else {
      if (comm->getRank() == 0) {
        std::cout << "Actual iters("<<actual_iters
              <<") != expected number of iterations ("
            <<expected_iters<<"), or resid-norm(" << norms[0] << ") >= 1.e-7"<<std::endl;
      }
    }
  }

  fpu_fix_end(&old_cw);

  timer.stop();
  if (comm->getRank() == 0) {
    std::cout << "proc 0 total program time: " << timer.totalElapsedTime()
        << std::endl;
  }

  return 0;
}

void process_command_line(int argc, char*argv[], std::string& xml_file)
{
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("xml_file", &xml_file, "XML Parameters file");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    throw std::runtime_error("Error parsing command-line.");
  }
}

#else
//else HAVE_IFPACK2_QD is not defined:

#include <iostream>
int main(int argc, char*argv[])
{
  std::cout << "belos_extprec_solve.cpp can't be compiled or run without the QD library." << std::endl;
  return 0;
}

#endif

