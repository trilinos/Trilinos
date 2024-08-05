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
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Stream.hpp"
#include "ROL_HelperFunctions.hpp"
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

    // Set algorithm parameters.
    ROL::ParameterList parlist;
    parlist.sublist("Step").set("Type", "Line Search");
    parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Newton-Krylov");
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);

    // Set initial guess.
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(2, 0.0);
    (*x_ptr)[0] = -1.2;
    (*x_ptr)[1] =  1.0;

    // Set scaling vector.
    ROL::Ptr<std::vector<RealT> > scale_ptr = ROL::makePtr<std::vector<RealT>>(2, 0.0);
    (*scale_ptr)[0] = 1.0;
    (*scale_ptr)[1] = 100.0;

    ROL::PrimalScaledStdVector<RealT> x(x_ptr, scale_ptr);

    // Define problem.
    ROL::OptimizationProblem<RealT> problem(ROL::makePtrFromRef(obj), ROL::makePtrFromRef(x));
    ROL::OptimizationSolver<RealT> solver(problem, parlist);
    // Solve problem.
    solver.solve(*outStream);

    // Set true solution.
    ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(2, 1.0);
    ROL::StdVector<RealT> xtrue(xtrue_ptr);

    // Compute dense Hessian.
    Teuchos::SerialDenseMatrix<int, RealT> H = ROL::computeDenseHessian(obj, x);
    Teuchos::SerialDenseMatrix<int, RealT> H_scaled = ROL::computeScaledDenseHessian(obj, x);
    *outStream << printMat(H);
    *outStream << printMat(H_scaled);

    // Compute Hessian error.
    Teuchos::SerialDenseMatrix<int, RealT> H_true(2, 2);
    H_true(0,0) = 802.0; H_true(0,1) = -400.0; H_true(1,0) = -400.0; H_true(1,1) = 200.0;
    Teuchos::SerialDenseMatrix<int, RealT> H_scaled_true(2, 2);
    H_scaled_true(0,0) = 802.0; H_scaled_true(0,1) = -400.0; H_scaled_true(1,0) = -4.0; H_scaled_true(1,1) = 2.0;
    H -= H_true;
    H_scaled -= H_scaled_true; H.normInf();
    if ((H.normInf() > sqrt(ROL::ROL_EPSILON<RealT>())) || (H_scaled.normInf() > sqrt(ROL::ROL_EPSILON<RealT>()))) {
      errorFlag += 1;
    }

    // Compute solution error.
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
