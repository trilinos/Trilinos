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
    \brief Shows how to use Belos in for a Krylov-Newton method    
    \author Created by G. von Winckel
*/


#include "ROL_TestObjectives.hpp"
#include "ROL_BelosKrylov.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <random>
#include <cstdlib>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);

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

      Teuchos::RCP<ROL::Step<RealT>> step;
      Teuchos::RCP<ROL::Objective<RealT>> obj;
 
      Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp(new Teuchos::ParameterList());
      std::string paramfile = "parameters.xml";
      Teuchos::updateParametersFromXmlFile(paramfile,Teuchos::Ptr<Teuchos::ParameterList>(&*parlist));

      // Get example number to use
      auto objNumber = static_cast<ROL::ETestObjectives>(parlist->get("Zoo Example Number",0));
 
      auto x0_rcp  = Teuchos::rcp(new std::vector<RealT>());
      auto x_rcp  = Teuchos::rcp(new std::vector<RealT>());
 
      ROL::StdVector<RealT> x(x_rcp);
      ROL::StdVector<RealT> x0(x0_rcp);
      
      ROL::getTestObjectives(obj,x0,x,objNumber);
      x.set(x0);

      // Make a Belos-Krylov solver if specified
      if(parlist->get("Use Belos",false)) { 
          Teuchos::RCP<ROL::Krylov<RealT>> krylov = Teuchos::rcp(new ROL::BelosKrylov<RealT>(*parlist));   
          step = Teuchos::rcp(new ROL::LineSearchStep<RealT>(krylov,*parlist));  
      }
      else { // Otherwise use ROL's default
          step = Teuchos::rcp(new ROL::LineSearchStep<RealT>(*parlist));
      }

      // Define Status Test
      RealT gtol  = 1e-12;  // norm of gradient tolerance
      RealT stol  = 1e-14;  // norm of step tolerance
      int   maxit = 100;    // maximum number of iterations
      ROL::StatusTest<RealT> status(gtol, stol, maxit);    

      // Define Algorithm
      ROL::DefaultAlgorithm<RealT> algo(*step,status,false);

      // Run Algorithm
      std::vector<std::string> output = algo.run(x, *obj, false);
      for ( unsigned i = 0; i < output.size(); i++ ) {
          std::cout << output[i];
      }

      // Check that the gradient norm is sufficiently small
      auto dim = x_rcp->size();
      auto g_rcp = Teuchos::rcp(new std::vector<RealT>(dim,0.0));  
      ROL::StdVector<RealT> g(g_rcp);
      RealT tol=0;
      obj->gradient(g,x,tol);

      if(g.norm()>gtol) {
          errorFlag = 1;
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

