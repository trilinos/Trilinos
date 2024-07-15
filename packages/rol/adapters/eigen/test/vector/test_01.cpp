// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_EigenVector.hpp"

#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>


int main(int argc, char *argv[]) {

  
  

  using RealT = double; 
  using E3V = ROL::Eigen3Vector<RealT>;
  using EigenVector = Eigen::Matrix<RealT,Eigen::Dynamic,1>;
 
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

  auto errtol = ROL::ROL_THRESHOLD<RealT>();

  // *** Test body.

  try {

    int dim = 10;

    auto x_ptr = ROL::makePtr<EigenVector>(dim);   
    auto y_ptr = ROL::makePtr<EigenVector>(dim);   
    auto z_ptr = ROL::makePtr<EigenVector>(dim);   

    E3V x(x_ptr);
    E3V y(y_ptr);
    E3V z(z_ptr); 

    //std::cout << x.dimension() << std::endl;


    RealT left = -1e0, right = 1e0;

    ROL::RandomizeVector(x, left, right );
    ROL::RandomizeVector(y, left, right );
    ROL::RandomizeVector(z, left, right );

    // Standard tests.
    auto consistency = x.checkVector(y, z, true, *outStream);
    ROL::StdVector<RealT> checkvec(ROL::makePtrFromRef(consistency));
    if (checkvec.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag++;
    }

    // Basis tests.
    // set x to first basis vector
    auto zp = x.basis(0);
    auto znorm = zp->norm();
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

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

