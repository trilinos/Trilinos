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

/*! \file  test_01.cpp
    \brief Test line search.
*/

#define USE_HESSVEC 1

#include "ROL_TestObjectives.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

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

  // *** Test body.

  try {

    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif

    for ( ROL::ETestOptProblem prob = ROL::TESTOPTPROBLEM_HS1; prob < ROL::TESTOPTPROBLEM_LAST; prob++ ) { 
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

      // Get Objective Function
      Teuchos::RCP<ROL::Vector<RealT> > x0, z;
      Teuchos::RCP<ROL::Objective<RealT> > obj;
      Teuchos::RCP<ROL::BoundConstraint<RealT> > con;
      ROL::getTestObjectives<RealT>(obj,con,x0,z,prob);
      Teuchos::RCP<ROL::Vector<RealT> > x = x0->clone();;

      // Get Dimension of Problem
      int dim = x0->dimension();
      parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

      // Error Vector
      Teuchos::RCP<ROL::Vector<RealT> > e = x0->clone();;
      e->zero();

      //ROL::ETrustRegion tr = ROL::TRUSTREGION_CAUCHYPOINT; 
      //ROL::ETrustRegion tr = ROL::TRUSTREGION_DOGLEG; 
      //ROL::ETrustRegion tr = ROL::TRUSTREGION_DOUBLEDOGLEG; 
      ROL::ETrustRegion tr = ROL::TRUSTREGION_TRUNCATEDCG; 
      parlist->sublist("Step").sublist("Trust Region").set("Subproblem Solver", ROL::ETrustRegionToString(tr));
      *outStream << std::endl << std::endl << ROL::ETrustRegionToString(tr) << std::endl << std::endl;
      
      // Define Algorithm
      ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);

      // Run Algorithm
      x->set(*x0);
      algo.run(*x, *obj, *con, true, *outStream);

      // Compute Error
      e->set(*x);
      e->axpy(-1.0,*z);
      *outStream << std::endl << "Norm of Error: " << e->norm() << std::endl;
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

