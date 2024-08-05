// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test ConicModel
*/

#include "ROL_OptimizationSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Rosenbrock.hpp"
#include "ROL_ConicApproximationModel.hpp"

#include <iostream>


int main(int argc, char *argv[]) {

  using RealT = double;
  using namespace ROL;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = makeStreamPtr( std::cout, argc > 1 );

  int errorFlag  = 0;

  // *** Test body.

  try {

    ZOO::getRosenbrock<RealT> rosenbrock;

    auto x   = rosenbrock.getInitialGuess();
    auto a   = x->clone();
    auto s   = x->clone();
    auto d   = x->clone();
    auto obj = rosenbrock.getObjective();

    RandomizeVector( *s );
    RandomizeVector( *a );
    RandomizeVector( *d );

    ConicApproximationModel<RealT> conic( obj, x, s, a );
    conic.checkGradient( *x, *d, true, *outStream );

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

