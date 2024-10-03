// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize the Zakharov function using NCG
*/

#define USE_HESSVEC 1

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_ParameterListConverters.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef std::vector<RealT>          vector;  
  typedef ROL::Vector<RealT>          V;      // Abstract vector
  typedef ROL::StdVector<RealT>       SV;     // Concrete vector containing std::vector data

  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  int errorFlag  = 0;

  // *** Example body.
 
  try {

    int dim = 10; // Set problem dimension. 

    std::string paramfile = "parameters.xml";
    auto parlist = ROL::getParametersFromXmlFile( paramfile );

   // Define algorithm.
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Iteration vector.
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim, 0.0);

    // Vector of natural numbers.
    ROL::Ptr<vector> k_ptr = ROL::makePtr<vector>(dim, 0.0);

    // For gradient and Hessian checks. 
    ROL::Ptr<vector> xtest_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> d_ptr     = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> v_ptr     = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> hv_ptr    = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> ihhv_ptr  = ROL::makePtr<vector>(dim, 0.0);
  

    RealT left = -1e0, right = 1e0; 
    for (int i=0; i<dim; i++) {
      (*x_ptr)[i]   = 2;
      (*k_ptr)[i]   = i+1.0;
    }

    ROL::Ptr<V> k = ROL::makePtr<SV>(k_ptr);
    SV x(x_ptr);

    // Check gradient and Hessian.
    SV xtest(xtest_ptr);
    SV d(d_ptr);
    SV v(v_ptr);
    SV hv(hv_ptr);
    SV ihhv(ihhv_ptr);

    ROL::RandomizeVector( xtest, left, right );
    ROL::RandomizeVector( d, left, right );
    ROL::RandomizeVector( v, left, right );

    ROL::ZOO::Objective_Zakharov<RealT> obj(k);

    obj.checkGradient(xtest, d, true, *outStream);                             *outStream << "\n"; 
    obj.checkHessVec(xtest, v, true, *outStream);                              *outStream << "\n";
    obj.checkHessSym(xtest, d, v, true, *outStream);                           *outStream << "\n";
   
    // Check inverse Hessian.
    RealT tol=0;
    obj.hessVec(hv,v,xtest,tol);
    obj.invHessVec(ihhv,hv,xtest,tol);
    ihhv.axpy(-1,v);
    *outStream << "Checking inverse Hessian" << std::endl;
    *outStream << "||H^{-1}Hv-v|| = " << ihhv.norm() << std::endl;
     

    // Run algorithm.
    algo.run(x, obj, true, *outStream);

    // Get True Solution
    ROL::Ptr<vector> xtrue_ptr = ROL::makePtr<vector>(dim, 0.0);
    SV xtrue(xtrue_ptr);

        
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr << std::endl;
    if ( abserr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
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

