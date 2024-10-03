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

#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_ObjectiveFromConstraint.hpp"

#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"


#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

   

  using Obj     = ROL::Objective<RealT>;
  using Con     = ROL::Constraint<RealT>;
  using V       = ROL::Vector<RealT>;
  using StdV    = ROL::StdVector<RealT>;
  using ScalarV = ROL::SingletonVector<RealT>;

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

  RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());

  // *** Test body.

  try {

    int xdim = 5;
    int cdim = 3;
    
    // Optimization vector
    ROL::Ptr<V> x  = ROL::makePtr<StdV>( ROL::makePtr<std::vector<RealT>>(xdim) );
    ROL::Ptr<V> c  = ROL::makePtr<StdV>( ROL::makePtr<std::vector<RealT>>(cdim));
    ROL::Ptr<V> e0 = c->basis(0);
    ROL::Ptr<V> e1 = c->basis(1);
    ROL::Ptr<V> e2 = c->basis(2);
 
    // Exact solution
    ROL::Ptr<V> sol   = x->clone();   
    ROL::Ptr<V> error = x->clone();

    ROL::Ptr<Obj> obj = ROL::nullPtr; 
    ROL::Ptr<Con> con = ROL::nullPtr;
    
    ROL::ZOO::getSimpleEqConstrained<RealT> SEC;
    obj = SEC.getObjective();
    con = SEC.getEqualityConstraint();
    x   = SEC.getInitialGuess();
    sol = SEC.getSolution();

    error->set(*sol);

    // Extract constraint components to make objectives
    ROL::Ptr<Obj> obj_0 = ROL::makePtr<ROL::ObjectiveFromConstraint<RealT>>( con, *e0 );
    ROL::Ptr<Obj> obj_1 = ROL::makePtr<ROL::ObjectiveFromConstraint<RealT>>( con, *e1 );
    ROL::Ptr<Obj> obj_2 = ROL::makePtr<ROL::ObjectiveFromConstraint<RealT>>( con, *e2 );

    // Create separate constraints from the objectives
    ROL::Ptr<Con> con_0 = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>( obj_0 );
    ROL::Ptr<Con> con_1 = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>( obj_1 );
    ROL::Ptr<Con> con_2 = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>( obj_2 );
    
    std::vector<ROL::Ptr<Con>> con_array;
    con_array.push_back(con_0);
    con_array.push_back(con_1);
    con_array.push_back(con_2);

    // Lagrange multipliers
    ROL::Ptr<V> l0 = ROL::makePtr<ScalarV>(0);
    ROL::Ptr<V> l1 = ROL::makePtr<ScalarV>(0);
    ROL::Ptr<V> l2 = ROL::makePtr<ScalarV>(0);
  
    std::vector<ROL::Ptr<V>> l_array;
    l_array.push_back(l0);
    l_array.push_back(l1);
    l_array.push_back(l2);
   
    ROL::OptimizationProblem<RealT> opt( obj,             // Objective
                                         x,               // Optimization vector
                                         ROL::nullPtr,   // No bound constraint
                                         con_array,       // Array of scalar equality constraints
                                         l_array);        // Array of scalar lagrange multipliers
 
    opt.check(*outStream);

    // Define algorithm.
    ROL::ParameterList parlist;
    std::string stepname = "Composite Step";
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);

    ROL::OptimizationSolver<RealT> solver( opt, parlist );

    solver.solve( *outStream );

    error->axpy(-1.0,*x);
    RealT error_norm = error->norm();
    
    *outStream << "\n\n Relative norm of final optimization vector error: " << error_norm << std::endl;
    
    if(error_norm>tol) 
      ++errorFlag;

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

