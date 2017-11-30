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

/*! \file  example_02.cpp
    \brief Test the OptimizationProblem interface on Step by taking
           alternating Linesearch and Trust-Region steps on the 
           Zahkarov function
*/

#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef std::vector<RealT>          vector;  
  typedef ROL::Vector<RealT>          V;      // Abstract vector
  typedef ROL::StdVector<RealT>       SV;     // Concrete vector containing std::vector data

  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.
 
  try {

    int dim = 10; // Set problem dimension. 

    ROL::Ptr<ParameterList> parlist = ROL::makePtr<ParameterList>();
    std::string paramfile = "parameters.xml";
    updateParametersFromXmlFile(paramfile,parlist.ptr());

    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim, 1.0);
    ROL::Ptr<vector> k_ptr = ROL::makePtr<vector>(dim, 0.0);

    ROL::Ptr<V> x = ROL::makePtr<SV>(x_ptr);  // Optimization vector
    ROL::Ptr<V> k = ROL::makePtr<SV>(k_ptr);  // Vector appearing in Zakharov objective

    ROL::Ptr<V> s = x->clone();            // Step vector

    for( int i=0; i<dim; ++i ) {
      (*k_ptr)[i] = i+1.0;
    }
    
    ROL::Ptr<ROL::Objective<RealT> > obj = ROL::makePtr<ROL::ZOO::Objective_Zakharov<RealT>>(k);
    
    ROL::OptimizationProblem<RealT> opt(obj,x);
    ROL::AlgorithmState<RealT> state;
   
    // Allocate iterate vector in algorithm state 
    state.iterateVec = x->clone();
    state.iterateVec->set(*x);
    state.minIterVec = x->clone();
 
    ROL::LineSearchStep<RealT>  ls(*parlist);
    ROL::TrustRegionStep<RealT> tr(*parlist);

    ls.initialize( opt, state );
    tr.initialize( opt, state );

    for( int iter = 0; iter<10; ++iter ) {
      ls.compute( *s, opt, state );
      ls.update( opt, *s, state );

      state.minIterVec->set(*x);
      state.minIter = state.iter;
      state.minValue = state.value;      

      *outStream << "LS fval = " << state.minValue << std::endl; 

      tr.compute( *s, opt, state );
      tr.update( opt, *s, state );

      state.minIterVec->set(*x);
      state.minIter = state.iter;
      state.minValue = state.value;      

      *outStream << "TR fval = " << state.minValue << std::endl; 
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



