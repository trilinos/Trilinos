// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve a finite minimax problem.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_BundleStep.hpp"
#include "ROL_BundleStatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Minimax2.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    // Initialize objective function.
    int dim = 2;
    ROL::ZOO::Minimax2<RealT> obj;

    // Initialize iteration vectors.
    ROL::Ptr<std::vector<RealT>> x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    (*x_ptr)[0] = 1.0; (*x_ptr)[1] = -0.1;
    ROL::StdVector<RealT> x(x_ptr);
    ROL::Ptr<std::vector<RealT>> z_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    (*z_ptr)[0] = 1.0; (*z_ptr)[1] = 1.0;
    ROL::StdVector<RealT> z(z_ptr);

    // Algorithmic input parameters.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::BundleStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::BundleStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Run algorithm.
    algo.run(x, obj, true, *outStream);

    // Compute error 
    ROL::Ptr<ROL::Vector<RealT>> diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,z);
    RealT error = diff->norm();
    *outStream << "\nAbsolute Error: " << error << "\n";
    *outStream <<   "Relative Error: " << error/z.norm() << "\n";
    errorFlag = ((error > 1e4*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);
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

