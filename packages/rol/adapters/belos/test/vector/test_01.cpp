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


#include "ROL_Zakharov.hpp"
#include "ROL_BelosKrylov.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <cstdlib>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    
      int dim = 10;

      ROL::Ptr<ROL::Step<RealT> > step;
 

      Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
      std::string paramfile = "parameters.xml";
      Teuchos::updateParametersFromXmlFile(paramfile,parlist.ptr());

      // Iteration Vector 
      ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);

      // Vector of natural numbers
      ROL::Ptr<std::vector<RealT> > k_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);

      for (int i=0; i<dim; i++) {
          (*x_ptr)[i]   = 4;
          (*k_ptr)[i]   = i+1.0;
       }

       ROL::Ptr<ROL::Vector<RealT> > k = ROL::makePtr<ROL::StdVector<RealT>>(k_ptr);
       ROL::StdVector<RealT> x(x_ptr);

       ROL::ZOO::Objective_Zakharov<RealT> obj(k);

      // Make a Belos-Krylov solver if specified
      if(parlist->get("Use Belos",false)) { 
          ROL::Ptr<ROL::Krylov<RealT> > krylov = ROL::makePtr<ROL::BelosKrylov<RealT>>(*parlist);   
          step = ROL::makePtr<ROL::LineSearchStep<RealT>>(*parlist,ROL::nullPtr,ROL::nullPtr,krylov);  
      }
      else { // Otherwise use ROL's default
          step = ROL::makePtr<ROL::LineSearchStep<RealT>>(*parlist);
      }

      // Define Status Test
      RealT gtol  = 1e-12;  // norm of gradient tolerance
      RealT stol  = 1e-14;  // norm of step tolerance
      int   maxit = 100;    // maximum number of iterations
      ROL::Ptr<ROL::StatusTest<RealT> > status = ROL::makePtr<ROL::StatusTest<RealT>>(gtol, stol, maxit);    

      // Define Algorithm
      ROL::Algorithm<RealT> algo(step,status,false);

      // Run Algorithm
      algo.run(x, obj, true, *outStream);

      // Get True Solution
      ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
      ROL::StdVector<RealT> xtrue(xtrue_ptr);
        
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

