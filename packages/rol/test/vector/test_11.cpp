// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_11.cpp
    \brief Test std::array interface.
*/

#include "ROL_StdArray.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

using RealT = double;

constexpr auto dim = 100u;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  ROL::StdArray<RealT,dim>::initialize_pool();

  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0) outStream = ROL::makePtrFromRef(std::cout);
  else            outStream = ROL::makePtrFromRef(bhs);

  auto print_pool_count = [&outStream]() {
    *outStream << "Currently using " << ROL::StdArray<RealT, dim>::pool_count() 
               << " Vectors from the pool" << std::endl;
  };

  int errorFlag  = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  // *** Test body.
 
  print_pool_count(); 

  try {


    ROL::StdArray<RealT, dim> x, y, z;

    print_pool_count(); 

    x.randomize();
    y.randomize();
    z.randomize();

    // Standard tests.
    auto consistency = x.checkVector(y, z, true, *outStream);
    ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );

    if (checkvec.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>())) errorFlag++;
    
    // Basis tests.
    // set x to first basis vector
    auto zp = x.clone();
    zp = x.basis(0);
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

    print_pool_count(); 
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  print_pool_count(); 

  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";


  return 0;

}

