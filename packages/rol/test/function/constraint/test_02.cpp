// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test ScalarLinearConstraint
           
*/

#include "ROL_ScalarLinearConstraint.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  using V = ROL::Vector<RealT>;

   

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

    int xdim = 5;
    
    // Optimization vector
    ROL::Ptr<V> a  = ROL::makePtr<ROL::StdVector<RealT>>( ROL::makePtr<std::vector<RealT>>(xdim) );
    ROL::Ptr<V> c  = ROL::makePtr<ROL::SingletonVector<RealT>>( 0.0 );

    ROL::Ptr<V> x = a->clone();
    ROL::Ptr<V> d = x->clone();
    ROL::Ptr<V> v = c->clone();

    RealT b = 0.5;
 
    ROL::RandomizeVector(*a);
    ROL::RandomizeVector(*x);
    ROL::RandomizeVector(*d);
    ROL::RandomizeVector(*v);

    std::cout << "a = "; a->print(*outStream); std::cout << std::endl;
    std::cout << "x = "; x->print(*outStream); std::cout << std::endl;
    std::cout << "d = "; d->print(*outStream); std::cout << std::endl;
    std::cout << "v = "; v->print(*outStream); std::cout << std::endl;

    ROL::ScalarLinearConstraint<RealT> con( a, b );

    con.checkApplyJacobian(*x,*d,*c,true,*outStream);
    con.checkAdjointConsistencyJacobian(*v,*d,*x,true,*outStream);
    con.checkApplyAdjointHessian(*x,*v,*d,*x,true,*outStream);

    

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

