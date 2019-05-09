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
    \brief Example of how to supply ROL with parameters from a JSON file.
           Requires that 
           <a href="https://github.com/open-source-parsers/jsoncpp">json-cpp</a> be installed.

           To build this example, add the following to the cmake call in your trilinos build script
           -D TPL_ENABLE_JSONCPP=ON \
           -D JSONCPP_INCLUDE_DIRS:PATH=/usr/include/jsoncpp \
           -D JSONCPP_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu \
           -D JSONCPP_LIBRARY_NAMES:STRING="jsoncpp" \

           These example paths above are default for Ubuntu 64 bit if jsoncpp is installed using 
           sudo apt-get install libjsoncpp-dev  
    
    \author Created by Greg von Winckel

*/

#include "ROL_Algorithm.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"

#include "example_01.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

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

  // *** Example body.

  try {

    int dim = 10;

    ROL::ParameterList parlist;

    std::string jsonFileName("parameters.json");
    parlist.setName("Imported from " + jsonFileName);
    
    // Load json parameters into a ROL::ParameterList  
    ROL::JSON_Parameters(jsonFileName,parlist);
    std::string stepname = "Trust Region"; // can we obtain this from parlist?  or jsonFile?

    ROL::StepFactory<RealT> stepFactory;
    ROL::Ptr<ROL::Step<RealT> > step = stepFactory.getStep(stepname, parlist);

    // Define Status Test
    RealT gtol  = parlist.get("Gradient Tolerance",1e-12); 
    RealT stol  = parlist.get("Step Tolerance",1e-14);  
    int   maxit = parlist.get("Maximum Number of Iterations",100); 
    ROL::Ptr<ROL::StatusTest<RealT> > status = ROL::makePtr<ROL::StatusTest<RealT>>(gtol, stol, maxit);           

    ROL::Algorithm<RealT> algo(step,status,false);

    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 1.0);
    ROL::Ptr<std::vector<RealT> > k_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);

    ROL::StdVector<RealT> x(x_ptr);
    ROL::Ptr<ROL::Vector<RealT> > k = ROL::makePtr<ROL::StdVector<RealT>>(k_ptr);

    for(int i=0;i<dim;++i) {
        (*k_ptr)[i] = i+1.0;
    }
 
    ROL::ZOO::Objective_Zakharov<RealT> obj(k);

    // Run Algorithm
    algo.run(x, obj, true, *outStream);

    // Get True Solution
    ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::StdVector<RealT> xtrue(xtrue_ptr);

    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    if ( abserr > sqrt(ROL::ROL_EPSILON<RealT>())) ) {
      errorFlag += 1;
    }

    // Make an XML file containing the supplied parameters
    Teuchos::writeParameterListToXmlFile(parlist,"parameters.xml");
 
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

