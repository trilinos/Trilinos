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

/*! \file  test_04.cpp
    \brief Test bound constrained trust-region steps.
*/

#define USE_HESSVEC 1

#include "ROL_GetTestProblems.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


#include <iostream>
//#include <fenv.h>

typedef double RealT;

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

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

  // *** Test body.

  try {

    std::string filename = "input.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif
    parlist->sublist("Step").set("Type","Trust Region");

    for ( ROL::ETestOptProblem prob = ROL::TESTOPTPROBLEM_ROSENBROCK; prob < ROL::TESTOPTPROBLEM_LAST; prob++ ) { 
      // Get Objective Function
      ROL::Ptr<ROL::Vector<RealT> > x0;
      std::vector<ROL::Ptr<ROL::Vector<RealT> > > z;
      ROL::Ptr<ROL::OptimizationProblem<RealT> > problem;
      ROL::GetTestProblem<RealT>(problem,x0,z,prob);
      if (problem->getProblemType() == ROL::TYPE_B) {
        if ( prob == ROL::TESTOPTPROBLEM_HS2 || prob == ROL::TESTOPTPROBLEM_BVP ) {
          parlist->sublist("Step").sublist("Line Search").set("Initial Step Size",1.e-4);
          parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",-1.e1);
          parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e-4);
          parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
        }
        else if ( prob == ROL::TESTOPTPROBLEM_HS25 ) {
          parlist->sublist("Step").sublist("Line Search").set("Initial Step Size",1.0);
          parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",1.e3);
          parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e4);
          parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
        }
        else {
          parlist->sublist("Step").sublist("Line Search").set("Initial Step Size",1.0);
          parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",-1.e1);
          parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e4);
          parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
        }
        parlist->sublist("General").set("Scale for Epsilon Active Sets",1.0);
        if ( prob == ROL::TESTOPTPROBLEM_HS4 ) {
          parlist->sublist("General").set("Scale for Epsilon Active Sets",1.e-2);
        }
        *outStream << std::endl << std::endl << ROL:: ETestOptProblemToString(prob)  << std::endl << std::endl;

        // Get Dimension of Problem
        int dim = x0->dimension();
        parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

        // Error Vector
        ROL::Ptr<ROL::Vector<RealT> > e = x0->clone();
        e->zero();

        //ROL::ETrustRegion tr = ROL::TRUSTREGION_CAUCHYPOINT; 
        //ROL::ETrustRegion tr = ROL::TRUSTREGION_DOGLEG; 
        //ROL::ETrustRegion tr = ROL::TRUSTREGION_DOUBLEDOGLEG; 
        ROL::ETrustRegion tr = ROL::TRUSTREGION_TRUNCATEDCG; 
        //ROL::ETrustRegion tr = ROL::TRUSTREGION_LINMORE; 
        parlist->sublist("Step").sublist("Trust Region").set("Subproblem Solver", ROL::ETrustRegionToString(tr));
        *outStream << std::endl << std::endl << ROL::ETrustRegionToString(tr) << std::endl << std::endl;

        // Define Solver
        ROL::OptimizationSolver<RealT> solver(*problem,*parlist);

        // Run Solver
        solver.solve(*outStream);

        // Compute Error
        RealT err(0);
        for (int i = 0; i < static_cast<int>(z.size()); ++i) {
          e->set(*x0);
          e->axpy(-1.0,*z[i]);
          if (i == 0) {
            err = e->norm();
          }
          else {
            err = std::min(err,e->norm());
          }
        }
        *outStream << std::endl << "Norm of Error: " << err << std::endl;
      }
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}
