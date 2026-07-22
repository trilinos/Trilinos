// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_05.cpp
    \brief Validate new augmented Lagrangian algorithm
*/

#include "ROL_HS32.hpp"
#include "ROL_Problem.hpp"
#include "ROL_TypeG_AugmentedLagrangianAlgorithm2.hpp"

#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  ROL::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    ROL::ZOO::getHS32<RealT> HS32;
    ROL::Ptr<ROL::Objective<RealT>>       obj   = HS32.getObjective();
    ROL::Ptr<ROL::Vector<RealT>>          x     = HS32.getInitialGuess();
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd   = HS32.getBoundConstraint();
    ROL::Ptr<ROL::Constraint<RealT>>      econ  = HS32.getEqualityConstraint();
    ROL::Ptr<ROL::Vector<RealT>>          emul  = HS32.getEqualityMultiplier();
    ROL::Ptr<ROL::Constraint<RealT>>      icon  = HS32.getInequalityConstraint();
    ROL::Ptr<ROL::Vector<RealT>>          imul  = HS32.getInequalityMultiplier();
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd  = HS32.getSlackBoundConstraint();

    ROL::Ptr<ROL::Problem<RealT>> problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);

    problem->addBoundConstraint(bnd);
    problem->addConstraint("inequality constraint",icon,imul,ibnd);
    problem->addConstraint("equality constraint",econ,emul);
    problem->finalize(false,true,std::cout,false);
    problem->check(true,*outStream,x);

    ROL::ParameterList list;
    list.sublist("General").set("Output Level", 3);
    list.sublist("Status Test").set("Iteration Limit", 20);
    list.sublist("Status Test").set("Constraint Tolerance", 1.e-4);
    list.sublist("Status Test").set("Step Tolerance", -1.);
    list.sublist("Step").set("Type", "Augmented Lagrangian 2");
    list.sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Iteration Limit", 100);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Initial Penalty Parameter", false);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Initial Penalty Parameter", 1.e1);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Penalty Parameter Growth Factor", 5.0);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Initial Optimality Tolerance", 1e-6);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Initial Feasibility Tolerance", 1e-5);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Composite Step");
    list.sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Fix Tolerance", true);
    list.sublist("Step").sublist("Composite Step").set("Initial Radius", 1.e-1);
    list.sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Nominal Relative Tolerance", 1.e-14);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling", false);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Level of Hessian Approximation", 1);
    // list.sublist("Step").sublist("Augmented Lagrangian").sublist("bounds").set("Penalty Parameter Growth Factor",1.0);
    ROL::Solver<RealT> solver(problem,list);
    solver.solve(*outStream);

    std::cout << std::endl << "Solution: ";
    x->print(*outStream);
    std::cout << std::endl << std::endl;
  }

  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
