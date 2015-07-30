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

/*! \file  example_02.cpp
*/

#include "ROL_Zakharov.hpp"
#include "example_02.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;  

  typedef ROL::Objective<RealT>          Obj;
  typedef ROL::StdVector<RealT>          V;

  typedef RCP<Obj>                       PObj;
  typedef RCP<V>                         PV;
 
  typedef std::vector<RealT>             Vec;

  typedef UnconstrainedTest<RealT>       Test;

   
  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy)
  // command-line argument is provided.
  int iprint     = argc - 1;
  RCP<std::ostream> outStream;
  oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = rcp(&std::cout, false);
  else
    outStream = rcp(&bhs, false);

  try {

    int maxit = 100; // Maximum number of iterations

    size_t dim = 10; // Set problem dimension. 

    RealT initVal;

    if(argc>1) {
      initVal = static_cast<RealT>( atof(argv[1]) );
    }
    else {
      initVal = 2.0;
    }


    RCP<ParameterList> parlist = rcp(new ParameterList());
    std::string paramfile = "parameters.xml";
    updateParametersFromXmlFile(paramfile,Ptr<ParameterList>(&*parlist));

    std::string Option1 = parlist->get("Option 1","Nonlinear CG Type");
    std::string Option2 = parlist->get("Option 2","Linesearch Type");


    RCP<Vec> x0_rcp = rcp( new Vec(dim,0.0) );
    RCP<Vec> x_rcp  = rcp( new Vec(dim,0.0) );
    RCP<Vec> k_rcp  = rcp( new Vec(dim,0.0) );

    for(size_t i=0; i<dim; ++i) {
      (*x0_rcp)[i] = initVal;
      (*k_rcp)[i] = i+1.0;
    }
    
    RCP<V> k = rcp( new V(k_rcp) );    // Vector of natural numbers

    PV x0 = rcp( new V(x0_rcp) );     // Initial optimization vector for reset
    
    PObj obj = rcp( new ROL::ZOO::Objective_Zakharov<RealT>(k) );
    
    Test test(obj,x0,*parlist,maxit);  // Create a test

    *outStream << "Printing final gradient norm for a variety of methods" << std::endl;
 
    test.run(Option1,Option2);

  }  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
  }; // end try
}

