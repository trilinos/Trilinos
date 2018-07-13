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

#ifndef HS_PROBLEM_038_HPP
#define HS_PROBLEM_038_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_038 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a1 = x[1]-x[0]*x[0];
    ScalarT a2 = 1.0-x[1];
    ScalarT a3 = x[3]-x[2]*x[2];
    ScalarT a4 = 1.0-x[2];
    ScalarT a5 = x[1]-1.0;
    ScalarT a6 = x[3]-1.0;

    return 100*a1*a1 + a2*a2 + 90*a3*a3 + a4*a4 + 10.1*(a5*a5+a6*a6) + 19.8*a5*a6;
 }
};
} // namespace HS_038

template<class Real> 
class Problem_038 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;

public:

  Problem_038() : NP( dimension_x() ) {
    NP::setLower(0, -10.0);
    NP::setUpper(0,  10.0);
    NP::setLower(1, -10.0);
    NP::setUpper(1,  10.0);
    NP::setLower(2, -10.0);
    NP::setUpper(2,  10.0);
    NP::setLower(3, -10.0);
    NP::setUpper(3,  10.0);

 }

  int dimension_x()  { return 4; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_038::Obj>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-3.0,-1.0,-3.0,-1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(19192.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = { 1.0, 1.0, 1.0, 1.0 };

    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_038_HPP
