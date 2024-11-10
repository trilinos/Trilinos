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
#include "ROL_HS42.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_TypeE_StabilizedLCLAlgorithm.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::string filename = "input.xml";
  auto parlist = ROL::getParametersFromXmlFile( filename );

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
    *outStream << std::endl << "Simple Equality Constrained Problem" << std::endl << std::endl;
    ROL::Ptr<ROL::Vector<RealT>> x;
    std::vector<ROL::Ptr<ROL::Vector<RealT>>> sol;
    ROL::Ptr<ROL::OptimizationProblem<RealT>> optProb;
    ROL::ZOO::getSimpleEqConstrained<RealT> SEC;
    SEC.get(optProb,x,sol);
    ROL::Ptr<ROL::Vector<RealT>> error = x->clone();

    // Solve optimization problem
    parlist->sublist("General").set("Output Level",iprint);
    parlist->sublist("Step").sublist("Stabilized LCL").set("Subproblem Iteration Limit",1000);
    parlist->sublist("Step").sublist("Stabilized LCL").set("Use Default Problem Scaling",false);
    parlist->sublist("Step").sublist("Stabilized LCL").set("Use Default Initial Penalty Parameter",false);
    parlist->sublist("Step").sublist("Stabilized LCL").set("Initial Penalty Parameter",10.0);
    parlist->sublist("Step").sublist("Stabilized LCL").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Status Test").set("Gradient Tolerance",1e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1e-12);
    parlist->sublist("Status Test").set("Constraint Tolerance",1e-8);
    ROL::TypeE::StabilizedLCLAlgorithm<RealT> algo(*parlist);
    algo.run(*optProb->getSolutionVector(),*optProb->getObjective(),
             *optProb->getConstraint(),*optProb->getMultiplierVector(),
             *outStream);

    error->set(*sol[0]);
    error->axpy(static_cast<RealT>(-1), *optProb->getSolutionVector());
    RealT solnErr = error->norm();

    *outStream << "Distance from true solution: " << solnErr << "\n";

    errorFlag += (solnErr < static_cast<RealT>(1e-6)) ? 0 : 1;

    *outStream << std::endl << "Hock and Schittkowski Problem #42" << std::endl << std::endl;
    ROL::ZOO::getHS42<RealT> HS42;
    ROL::Ptr<ROL::Vector<RealT>>     xvec       = HS42.getInitialGuess();
    ROL::Ptr<ROL::Objective<RealT>>  obj        = HS42.getObjective();
    ROL::Ptr<ROL::Constraint<RealT>> linear_con = ROL::makePtr<ROL::ZOO::Constraint_HS42a<RealT>>();
    ROL::Ptr<ROL::Vector<RealT>>     linear_mul = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    ROL::Ptr<ROL::Constraint<RealT>> con        = ROL::makePtr<ROL::ZOO::Constraint_HS42b<RealT>>();
    ROL::Ptr<ROL::Vector<RealT>>     mul        = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    ROL::TypeE::StabilizedLCLAlgorithm<RealT> algo1(*parlist);
    algo1.run(*xvec,*obj,*con,*mul,*linear_con,*linear_mul,*outStream);

    error = HS42.getSolution();
    error->axpy(static_cast<RealT>(-1), *xvec);
    RealT solnErr1 = error->norm();

    *outStream << "Distance from true solution: " << solnErr1 << "\n";

    errorFlag += (solnErr1 < static_cast<RealT>(1e-6)) ? 0 : 1;

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
