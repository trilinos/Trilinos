// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test stacktrace
*/

#include <exception>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_stacktrace.hpp"
#include "ROL_Stream.hpp"

int main( int argc, char* argv[] ) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  try {
    ROL_TEST_FOR_EXCEPTION( 0<1, std::logic_error, "End Result: TEST PASSED\n")
  }
  catch( std::exception& e ) {
    *outStream << e.what() << std::endl;   
  }

  return 0;
}

