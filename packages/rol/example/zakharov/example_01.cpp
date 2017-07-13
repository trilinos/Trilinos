// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize the Zakharov function using NCG
*/

#define USE_HESSVEC 1

#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"
#include "ROL_ParameterListConverters.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef std::vector<RealT>          vector;  
  typedef ROL::Vector<RealT>          V;      // Abstract vector
  typedef ROL::StdVector<RealT>       SV;     // Concrete vector containing std::vector data

  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.
 
  try {

    int dim = 10; // Set problem dimension. 

    RCP<ParameterList> parlist = rcp(new ParameterList());
    std::string paramfile = "parameters.xml";
    updateParametersFromXmlFile(paramfile,parlist.ptr());

   // Define algorithm.
    ROL::Algorithm<RealT> algo("Trust-Region",*parlist);

    // Iteration vector.
    RCP<vector> x_rcp = rcp( new vector(dim, 0.0) );

    // Vector of natural numbers.
    RCP<vector> k_rcp = rcp( new vector(dim, 0.0) );

    // For gradient and Hessian checks. 
    RCP<vector> xtest_rcp = rcp( new vector(dim, 0.0) );
    RCP<vector> d_rcp     = rcp( new vector(dim, 0.0) );
    RCP<vector> v_rcp     = rcp( new vector(dim, 0.0) );
    RCP<vector> hv_rcp    = rcp( new vector(dim, 0.0) );
    RCP<vector> ihhv_rcp  = rcp( new vector(dim, 0.0) );
  

    RealT left = -1e0, right = 1e0; 
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i]   = 2;
      (*k_rcp)[i]   = i+1.0;
    }

    RCP<V> k = rcp(new SV(k_rcp) );
    SV x(x_rcp);

    // Check gradient and Hessian.
    SV xtest(xtest_rcp);
    SV d(d_rcp);
    SV v(v_rcp);
    SV hv(hv_rcp);
    SV ihhv(ihhv_rcp);

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
    RCP<vector> xtrue_rcp = rcp( new vector(dim, 0.0) );
    SV xtrue(xtrue_rcp);

        
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr << std::endl;
    if ( abserr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

