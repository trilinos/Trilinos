// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test scalar minimization algorithms on test 01.
*/
#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ScalarMinimizationTest01.hpp"

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

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // *** Test body.

  try {

    ROL::ParameterList parlist;
    parlist.sublist("Scalar Minimization").set("Type","Brent's");
    parlist.sublist("Scalar Minimization").sublist("Brent's").set("Tolerance",1.e-10);
    parlist.sublist("Scalar Minimization").sublist("Brent's").set("Iteration Limit",1000);
    ROL::Ptr<ROL::ScalarMinimizationTest<RealT> > smt
      = ROL::makePtr<ROL::ScalarMinimizationTest01<RealT>>(parlist);
    *outStream << "\nBrent's Method\n";
    bool flag = smt->test(*outStream);
    errorFlag += (flag ? 0 : 1);

    parlist.sublist("Scalar Minimization").set("Type","Bisection");
    parlist.sublist("Scalar Minimization").sublist("Bisection").set("Tolerance",1.e-10);
    parlist.sublist("Scalar Minimization").sublist("Bisection").set("Iteration Limit",1000);
    smt = ROL::makePtr<ROL::ScalarMinimizationTest01<RealT>>(parlist);
    *outStream << "\nBisection Method\n";
    flag = smt->test(*outStream);
    errorFlag += (flag ? 0 : 1);

    parlist.sublist("Scalar Minimization").set("Type","Golden Section");
    parlist.sublist("Scalar Minimization").sublist("Golden Section").set("Tolerance",1.e-10);
    parlist.sublist("Scalar Minimization").sublist("Golden Section").set("Iteration Limit",1000);
    smt = ROL::makePtr<ROL::ScalarMinimizationTest01<RealT>>(parlist);
    *outStream << "\nGolden Section Method\n";
    flag = smt->test(*outStream);
    errorFlag += (flag ? 0 : 1);
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

