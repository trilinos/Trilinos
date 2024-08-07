// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_16.cpp
    \brief Test augmented Lagrangian step.
*/

#include "ROL_GetTestProblems.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

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
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",10.0);

    // Setup optimization problem
    ROL::Ptr<ROL::Vector<RealT>> x0;
    std::vector<ROL::Ptr<ROL::Vector<RealT>>> z;
    ROL::Ptr<ROL::ZOO::getCubic<RealT>> cubic;
    ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    for (int i = 0; i < 3; ++i) {
      cubic = ROL::makePtr<ROL::ZOO::getCubic<RealT>>(i);
      cubic->get(problem,x0,z);

      // Check Derivatives
      problem->check(*outStream);

      // Setup optimization solver
      ROL::OptimizationSolver<RealT> solver(*problem,*parlist);
      solver.solve(*outStream);

      // Compute Error
      ROL::Ptr<ROL::Vector<RealT>> e = x0->clone();
      e->zero();
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
      RealT tol = static_cast<RealT>(1e-3)*std::max(z[0]->norm(),static_cast<RealT>(1));
      errorFlag += ((err < tol) ? 0 : 1);
    }
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

