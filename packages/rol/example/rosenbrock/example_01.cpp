// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize Rosenbrock's function using Newton-Krylov.
    \addtogroup examples_group
*/

#define USE_HESSVEC 1

#include "ROL_Rosenbrock.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Stream.hpp"
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

  // *** Example body.

  try {

    ROL::ZOO::Objective_Rosenbrock<RealT> obj;
    int dim = 100; // Set problem dimension. Must be even.

    // Set parameters.
    ROL::ParameterList parlist;
    parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Newton-Krylov");
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);

    // Define algorithm.
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::LineSearchStep<RealT>>(parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Iteration Vector
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    // Set Initial Guess
    for (int i=0; i<dim/2; i++) {
      (*x_ptr)[2*i]   = -1.2;
      (*x_ptr)[2*i+1] =  1.0;
    }
    ROL::StdVector<RealT> x(x_ptr);

    // Run Algorithm
    algo.run(x, obj, true, *outStream);

    // Get True Solution
    ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(dim, 1.0);
    ROL::StdVector<RealT> xtrue(xtrue_ptr);
    
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    RealT relerr = abserr/xtrue.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
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

