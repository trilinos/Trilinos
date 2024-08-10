// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test the NonlinearProblem interface with the Hock & Schittkowski
           problem set
*/

#define OPTIMIZATION_PROBLEM_REFACTOR 

#include "HS_ProblemFactory.hpp"

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_ValidParameters.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iomanip>

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace ROL;
  
  
  
  typedef Teuchos::ParameterList PL;

  typedef Vector<RealT>              V;
  typedef OptimizationProblem<RealT> OPT;
  typedef NonlinearProgram<RealT>    NLP;
  typedef AlgorithmState<RealT>      STATE;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag   = 0;
  int numProblems = 40;

  RealT errtol = 1e-5; // Require norm of computed solution error be less than this

  std::vector<int> passedTests;
  std::vector<int> failedTests;
  std::vector<std::string> converged;
  std::vector<std::string> stepname;

  HS::ProblemFactory<RealT> factory;

  ROL::Ptr<NLP>  nlp;
  ROL::Ptr<OPT>  opt;

  // Get two copies so we can validate without setting defaults
  ROL::Ptr<PL>   parlist = ROL::makePtr<PL>();
  ROL::Ptr<PL>   parlistCopy = ROL::makePtr<PL>();
  ROL::Ptr<PL>   test = ROL::makePtr<PL>();

  Teuchos::updateParametersFromXmlFile(std::string("hs_parameters.xml"),parlist.ptr());
  Teuchos::updateParametersFromXmlFile(std::string("hs_parameters.xml"),parlistCopy.ptr());
  Teuchos::updateParametersFromXmlFile(std::string("test_parameters.xml"),test.ptr());

  ROL::Ptr<const PL> validParameters = getValidROLParameters();

  parlistCopy->validateParametersAndSetDefaults(*validParameters);

  ROL::Ptr<V>           x;
  ROL::Ptr<const STATE> algo_state;

  bool problemSolved;

  // *** Test body.
  try {
 
    for(int n=1; n<=numProblems; ++n) {
      *outStream << "\nHock & Schittkowski Problem " << std::to_string(n) << std::endl;

      nlp = factory.getProblem(n);
      opt = nlp->getOptimizationProblem();
      EProblem problemType = opt->getProblemType();

      std::string str;

      switch( problemType ) {

        case TYPE_U:                                   
          str = test->get("Type-U Step","Trust Region");
        break;

        case TYPE_B:           
          str = test->get("Type-B Step","Trust Region");
        break;

        case TYPE_E:  
          str = test->get("Type-E Step","Composite Step");
        break;

        case TYPE_EB:
          str = test->get("Type-EB Step","Augmented Lagrangian");
        break;

        case TYPE_LAST:
          ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Error: Unsupported problem type!");
        break;
      }

      stepname.push_back(str);

      parlist->sublist("Step").set("Type",str);

      OptimizationSolver<RealT> solver( *opt, *parlist );

      solver.solve(*outStream);

      algo_state = solver.getAlgorithmState();

      x = opt->getSolutionVector();       
     
      problemSolved = nlp->foundAcceptableSolution( *x, errtol );

      RealT tol = std::sqrt(ROL_EPSILON<RealT>()); 

      *outStream << "Target objective value   = " << nlp->getSolutionObjectiveValue()   << std::endl;
      *outStream << "Attained objective value = " << opt->getObjective()->value(*x,tol) << std::endl;

      if( problemSolved ) {
        *outStream << "Computed an acceptable solution" << std::endl;        
        passedTests.push_back(n);
        converged.push_back(std::to_string(algo_state->nfval));
      }
      else {
        *outStream << "Failed to converge" << std::endl;
        failedTests.push_back(n);
        converged.push_back("Failed");
      }

    } // end for loop over HS problems  


    *outStream << "\nTests passed: ";
    for( auto const& value : passedTests ) {
      *outStream << value << " ";
    }
    *outStream << "\nTests failed: ";
    for( auto const& value : failedTests ) {
      *outStream << value << " ";
    }

    *outStream << std::endl;

    if( passedTests.size() > failedTests.size() ) {
      *outStream << "Most tests passed." << std::endl; 
    }
    else {
      *outStream << "Most tests failed." << std::endl;
      errorFlag++;
    }
         
    *outStream << "\n\nPERFORMANCE SUMMARY:\n\n";
    *outStream << std::setw(16) << "H&S Problem #" 
               << std::setw(24) << "Step Used" 
               << std::setw(12) << "#fval" << std::endl;
    *outStream << std::string(52,'-') << std::endl;

    for(int n=0; n<numProblems; ++n) {
      *outStream << std::setw(16) << n+1
                 << std::setw(24) << stepname[n] 
                 << std::setw(12) << converged[n] << std::endl;
    }
   

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
} 


