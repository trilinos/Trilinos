// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test std::vector interface.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;
typedef double ElementT;

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

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  // *** Test body.

  try {

    // Number of elements in StdVector objects
    int dim = 100; 

    ROL::StdVector<RealT, ElementT> x(dim);
    ROL::StdVector<RealT, ElementT> y(dim);
    ROL::StdVector<RealT, ElementT> z(dim);

    // Lower and upper bounds on randomized vector element values
    RealT left = -1e0, right = 1e0;

    // set x,y,z
    x.randomize(left,right);
    y.randomize(left,right);
    z.randomize(left,right);

    // Standard tests.
    auto consistency = x.checkVector(y, z, true, *outStream);

    ROL::StdVector<RealT, ElementT> checkvec( ROL::makePtrFromRef(consistency) );

    if (checkvec.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>())) errorFlag++;

    // Basis tests.
    // set x to first basis vector
    auto zp = x.basis(0);
    RealT znorm = zp->norm();
    *outStream << "Norm of ROL::Vector z (first basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // set x to middle basis vector
    zp = x.basis(dim/2);
    znorm = zp->norm();
    *outStream << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";

    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // set x to last basis vector
    zp = x.basis(dim-1);
    znorm = zp->norm();
    *outStream << "\nNorm of ROL::Vector z (last basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // Repeat the checkVector tests with a zero vector.
    x.scale(0.0);
    consistency = x.checkVector(x, x, true, *outStream);
    if (checkvec.norm() > 0.0) {
      errorFlag++;
    }
  }

  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;

}

