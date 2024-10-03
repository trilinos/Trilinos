// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize Rosenbrock's function.
    \addtogroup examples_group
*/

#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "example_01.hpp"

#include <iostream>

typedef double RealT;
typedef float  ElementT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFRomRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {

    Objective_Rosenbrock_AF<RealT, ElementT> obj;
    int dim = 10000; // Set problem dimension. Must be even.

    // Set parameters.
    std::string filename;
    if (sizeof(RealT) == sizeof(double)) {
      filename = "input-double.xml";
    }
    if (sizeof(RealT) == sizeof(float)) {
      filename = "input-float.xml";
    }
    *outStream << std::endl << "Using input file: " << filename << std::endl << std::endl;
    auto parlist = ROL::makePtr<ROL::ParameterList>()
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Define algorithm.
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::LineSearchStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Iteration Vector
    /***** Display ArrayFire info. *****/
    af::dtype afType = f32;
    if (sizeof(ElementT) == sizeof(double)) {
      afType = f64;
    }
    int device = 0;
    af::setDevice(device);
    std::string afinfo;
    afinfo = af::infoString(true);
    *outStream << std::endl << afinfo << std::endl;

    /***** Define ROL::Ptr to AF array with initial guess. *****/
    ElementT zero(0);
    ElementT op2(1.2);
    ElementT one(1);
    ROL::Ptr<af::array> x_ptr = ROL::makePtr<af::array>(dim,afType);
    *x_ptr = af::constant(zero,x_ptr->dims(),afType);
    // Set Initial Guess
    for (dim_t i=0; i<dim/2; i++) {
      (*x_ptr)(2*i)   = -op2;
      (*x_ptr)(2*i+1) =  one;
    }
    //af_print(*x_ptr);
    ROL::ArrayFireVector<RealT,ElementT> x(x_ptr);

    // Run Algorithm
    algo.run(x, obj, true, *outStream);
    //af_print(*x_ptr);

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

