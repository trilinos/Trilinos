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

#ifndef HS_PROBLEM_025_HPP
#define HS_PROBLEM_025_HPP

#include "ROL_NonlinearProgram.hpp"


namespace HS {

namespace HS_025 {
template<class Real>
class Obj {
private:
  template<class ScalarT>
  ScalarT f(int i, const std::vector<ScalarT> &x ) {
    Real u = 25.0 + std::pow(-50.0*std::log(0.01*i),2.0/3.0);
    return -0.01*i + std::exp( -std::pow(u-x[1],x[2]) /x[0] );
  }
public:

  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a=0;
    for(int i=1;i<100;++i) {
      ScalarT b = f(i,x);
      a  = a+b*b;
    }
    return a;
  }
};
}

template<class Real> 
class Problem_025 : public ROL::NonlinearProgram<Real> {
 
  

  typedef ROL::NonlinearProgram<Real>    NP;
  typedef ROL::Vector<Real>              V;
  typedef ROL::Objective<Real>           OBJ;

private:
public:

  Problem_025() : NP( dimension_x() ) {
    NP::setLower(0,0.1);
    NP::setUpper(0,100.0);

    NP::setLower(1,0.0);
    NP::setUpper(1,25.6);

    NP::setLower(2,0.0);
    NP::setUpper(2,5.0);
  }

  int dimension_x() { return 3; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_025::Obj>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {100,12.5,1.5};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(32.835);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    const Real x[] = {50.0,25.0,1.5};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

}

#endif // HS_PROBLEM_025_HPP
