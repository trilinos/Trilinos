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
#include "ROL_StatusTest.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef ROL::StdVector<RealT>       V;
  typedef std::vector<RealT>          Vec;

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
    updateParametersFromXmlFile(paramfile,Ptr<ParameterList>(&*parlist));

    parlist->set("Descent Type", "Newton-Krylov");

    // Define Step
    ROL::LineSearchStep<RealT> step(*parlist);

    // Define Status Test
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-13;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    ROL::StatusTest<RealT> status(gtol, stol, maxit);    

    // Define Algorithm
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);

    // Iteration Vector 
    RCP<Vec> x_rcp = rcp( new Vec(dim, 0.0) );

    // Vector of natural numbers
    RCP<Vec> k_rcp = rcp( new Vec(dim, 0.0) );

    // For gradient and Hessian checks 
    RCP<Vec> xtest_rcp = rcp( new Vec(dim, 0.0) );
    RCP<Vec> d_rcp     = rcp( new Vec(dim, 0.0) );
    RCP<Vec> v_rcp     = rcp( new Vec(dim, 0.0) );
    RCP<Vec> hv_rcp    = rcp( new Vec(dim, 0.0) );
    RCP<Vec> ihhv_rcp  = rcp( new Vec(dim, 0.0) );
  
    RealT left = -1e0, right = 1e0; 
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i]   = 2;
      (*k_rcp)[i]   = i+1.0;

      (*xtest_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*d_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*v_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    Teuchos::RCP<ROL::Vector<RealT> > k = Teuchos::rcp(new ROL::StdVector<RealT> (k_rcp) );
    ROL::StdVector<RealT> x(x_rcp);

    // Check gradient and Hessian
    V xtest(xtest_rcp);
    V d(d_rcp);
    V v(v_rcp);
    V hv(hv_rcp);
    V ihhv(ihhv_rcp);

    ROL::ZOO::Objective_Zakharov<RealT> obj(k);

    obj.checkGradient(xtest, d, true, *outStream);                             *outStream << "\n"; 
    obj.checkHessVec(xtest, v, true, *outStream);                              *outStream << "\n";
    obj.checkHessSym(xtest, d, v, true, *outStream);                           *outStream << "\n";
   
    // Check inverse Hessian 
    RealT tol=0;
    obj.hessVec(hv,v,xtest,tol);
    obj.invHessVec(ihhv,hv,xtest,tol);
    ihhv.axpy(-1,v);
    std::cout << "Checking inverse Hessian" << std::endl;
    std::cout << "||H^{-1}Hv-v|| = " << ihhv.norm() << std::endl;
     

    // Run Algorithm
    std::vector<std::string> output = algo.run(x, obj, false);
    for ( unsigned i = 0; i < output.size(); i++ ) {
      std::cout << output[i];
    }

    // Get True Solution
    RCP<Vec> xtrue_rcp = rcp( new Vec(dim, 0.0) );
    V xtrue(xtrue_rcp);

        
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    if ( abserr > sqrt(ROL::ROL_EPSILON) ) {
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

