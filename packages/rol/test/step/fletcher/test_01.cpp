// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test ConstraintFromObjective and ObjectiveFromConstraint
           which will help test ConstraintManager and OptimizationProblem
*/

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Fletcher.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::string filename = "input.xml";
  
  auto parlist = ROL::getParametersFromXmlFile( filename );

  using Opt     = ROL::OptimizationProblem<RealT>;
  using V       = ROL::Vector<RealT>;

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // *** Test body.

  try {

    // Set up optimization problem
    ROL::Ptr<V> x;
    std::vector<ROL::Ptr<V> > sol;
    ROL::Ptr<Opt> optProb;
    ROL::ZOO::getSimpleEqConstrained<RealT> SEC;
    SEC.get( optProb, x, sol );
    ROL::Ptr<V> error = x->clone();

    // Solve optimization problem
    ROL::OptimizationSolver<RealT> optSolver(*optProb, *parlist);
    optSolver.solve(*outStream);

    error->set(*sol[0]);
    error->axpy(static_cast<RealT>(-1), *x);
    RealT solnErr = error->norm();

    *outStream << "Distance from true solution: " << solnErr << "\n";

    errorFlag += (solnErr < static_cast<RealT>(1e-6)) ? 0 : 1;
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;
}
