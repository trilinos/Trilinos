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

#include "ROL_Zakharov.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

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

    Teuchos::ParameterList parlist;
    // Enumerations
    parlist.set("Descent Type",                           "Nonlinear-CG");
    parlist.set("Linesearch Type",                        "Cubic Interpolation");
    parlist.set("Linesearch Curvature Condition",         "Wolfe");
    // Linesearch Parameters
    parlist.set("Maximum Number of Function Evaluations", 20);
    parlist.set("Sufficient Decrease Parameter",          1.e-4);
    parlist.set("Curvature Conditions Parameter",         0.9);
    parlist.set("Backtracking Rate",                      0.5);
    parlist.set("Initial Linesearch Parameter",           1.0);
    parlist.set("User Defined Linesearch Parameter",      false);
    // Krylov Parameters
    parlist.set("Absolute Krylov Tolerance",              1.e-4);
    parlist.set("Relative Krylov Tolerance",              1.e-2);
    parlist.set("Maximum Number of Krylov Iterations",    10);
    // Define Step
    ROL::LineSearchStep<RealT> step(parlist);

    // Define Status Test
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-14;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    ROL::StatusTest<RealT> status(gtol, stol, maxit);    

    // Define Algorithm
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);

    // Iteration Vector 
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );

    // Vector of natural numbers
    Teuchos::RCP<std::vector<RealT> > k_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );

    // For gradient and Hessian checks 
    Teuchos::RCP<std::vector<RealT> > xtest_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > d_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > v_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > hv_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > ihhv_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
  
    RealT left = -1e0, right = 1e0; 
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i]   = 4;
      (*k_rcp)[i]   = i+1.0;

      (*xtest_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*d_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*v_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    Teuchos::RCP<ROL::Vector<RealT> > k = Teuchos::rcp(new ROL::StdVector<RealT> (k_rcp) );
    ROL::StdVector<RealT> x(x_rcp);

    // Check gradient and Hessian
    ROL::StdVector<RealT> xtest(xtest_rcp);
    ROL::StdVector<RealT> d(d_rcp);
    ROL::StdVector<RealT> v(v_rcp);
    ROL::StdVector<RealT> hv(hv_rcp);
    ROL::StdVector<RealT> ihhv(ihhv_rcp);

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
    Teuchos::RCP<std::vector<RealT> > xtrue_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> xtrue(xtrue_rcp);

        
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

