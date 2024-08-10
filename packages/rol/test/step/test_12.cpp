// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test line search.
*/

#define USE_HESSVEC 0

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Step.hpp"
#include "ROL_GetTestProblems.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;
  ROL::Ptr<ROL::Vector<RealT> > x0;
  std::vector<ROL::Ptr<ROL::Vector<RealT> > > x;
  ROL::Ptr<ROL::OptimizationProblem<RealT> > problem;
  ROL::GetTestProblem<RealT>(problem,x0,x,ROL::TESTOPTPROBLEM_HS1);
  ROL::Ptr<ROL::Objective<RealT> > obj = problem->getObjective();
  ROL::Ptr<ROL::BoundConstraint<RealT> > bnd = problem->getBoundConstraint();
  ROL::AlgorithmState<RealT> algo_state;

  // *** Test body.
  ROL::Step<RealT> step;
  int thrown = 0;
  try {
    try {
      step.compute(*x0,*x[0],*obj,*bnd,algo_state);
    }
    catch (ROL::Exception::NotImplemented& exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.update(*x0,*x[0],*obj,*bnd,algo_state);
    }
    catch (ROL::Exception::NotImplemented& exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.printName();
    }
    catch (ROL::Exception::NotImplemented& exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.printHeader();
    }
    catch (ROL::Exception::NotImplemented& exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.print(algo_state,true);
    }
    catch (ROL::Exception::NotImplemented& exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    errorFlag = (thrown==5) ? 0 : 1;
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}

