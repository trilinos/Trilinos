// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_10.cpp
    \brief Test bound constrainted interior point.
*/

#define USE_HESSVEC 1

#include "ROL_GetTestProblems.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Stream.hpp"
#include "ROL_InteriorPointStep.hpp"

#include "Teuchos_GlobalMPISession.hpp"

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

  // *** Test body.

  try {
    // Get Objective Function
    ROL::Ptr<ROL::Vector<RealT> > x0;
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > z;
    ROL::Ptr<ROL::OptimizationProblem<RealT> > problem;
    ROL::GetTestProblem<RealT>(problem,x0,z,ROL::TESTOPTPROBLEM_HS38);
    // Parse input
    std::string filename = "input.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif

    RealT mu = 0.1;            // Initial penalty parameter
    RealT factor = 0.1;        // Penalty reduction factor

    // Set solver parameters
    parlist->sublist("Step").set("Type","Interior Point");
    parlist->sublist("Step").sublist("Interior Point").set("Initial Barrier Penalty",mu);
    parlist->sublist("Step").sublist("Interior Point").set("Minimum Barrier Penalty",1e-8);
    parlist->sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor",factor);
    parlist->sublist("Step").sublist("Interior Point").set("Subproblem Iteration Limit",30);

    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist->sublist("Step").sublist("Composite Step").set("Output Level",0);

    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Iteration Limit",100);

    // Solve optimization problem with interior points
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::InteriorPointStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);
    algo.run(optProb, true, *outStream);

    // Compute Error
    ROL::Ptr<ROL::Vector<RealT> > e = x0->clone();
    RealT err(0);
    for (int i = 0; i < static_cast<int>(z.size()); ++i) {
      e->set(*x0);
      e->axpy(-1.0,*z[i]);
      if (i == 0) {
        err = e->norm();
      }
      else {
        err = std::min(err,e->norm());
      }
    }
    *outStream << std::endl << "Norm of Error: " << err << std::endl;
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

