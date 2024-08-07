// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test sketching QR.
*/


#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "ROL_Sketch.hpp"
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

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // *** Test body.
  try {

    int nrow = 128, ncol = 100, rank = 1, testrank = 6;
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(nrow, 0.0);
    ROL::StdVector<RealT> x(x_ptr);

    ROL::Sketch<RealT> sketch(x,ncol,rank,ROL::ROL_EPSILON<RealT>());
    bool flag = true;
    *outStream << std::endl;
    for (int i = 0; i < testrank; ++i) {
      *outStream << "Rank = " << i+1 << "  Test Rank = " << testrank << std::endl;
      sketch.setRank(i+1);
      flag = sketch.test(testrank,*outStream,1);
    }
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

