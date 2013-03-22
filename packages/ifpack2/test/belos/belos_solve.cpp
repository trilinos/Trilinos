/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#include "Ifpack2_Parameters.hpp"

#include "Tpetra_DefaultPlatform.hpp"

#include "build_problem.hpp"
#include "build_solver.hpp"


void process_command_line(int argc, char*argv[], std::string& xml_file);


int main(int argc, char*argv[])
{
  Teuchos::GlobalMPISession mpisess(&argc, &argv);

  bool success = true;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    Teuchos::Time timer("total");
    timer.start();

    Tpetra::DefaultPlatform::DefaultPlatformType& platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();

    typedef double Scalar;
    typedef int LO; //LocalOrdinal
    typedef int GO; //GlobalOrdinal
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
    typedef Tpetra::Operator<Scalar,LO,GO,Node>    TOP;
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

    *out << "Every proc reading parameters from xml_file: "
         << xml_file << std::endl;
    Teuchos::ParameterList test_params =
      Teuchos::ParameterXMLFileReader(xml_file).getParameters();

    //The build_problem function is located in build_problem.hpp.
    //Note that build_problem calls build_precond and sets a preconditioner on the
    //linear-problem, if a preconditioner is specified.

    Teuchos::RCP<BLinProb> problem =
      build_problem<Scalar,LO,GO,Node>(test_params, comm, platform.getNode());

    //The build_solver function is located in build_solver.hpp:

    Teuchos::RCP<BSolverMgr> solver = build_solver<Scalar,TMV,TOP>(test_params, problem);

    Belos::ReturnType ret = solver->solve();

    *out << "Converged in " << solver->getNumIters() << " iterations." << std::endl;

    Teuchos::RCP<TMV> R = Teuchos::rcp(new TMV(*problem->getRHS()));
    problem->computeCurrResVec(&*R, &*problem->getLHS(), &*problem->getRHS());
    Teuchos::Array<Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(R->getNumVectors());
    R->norm2(norms);

    if (norms.size() < 1) {
      throw std::runtime_error("ERROR: norms.size()==0 indicates R->getNumVectors()==0.");
    }

    *out << "2-Norm of 0th residual vec: " << norms[0] << std::endl;

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


void process_command_line(int argc, char*argv[], std::string& xml_file)
{
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("xml_file", &xml_file, "XML Parameters file");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    throw std::runtime_error("Error parsing command-line.");
  }
}

