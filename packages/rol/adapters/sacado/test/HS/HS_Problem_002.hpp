// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_002_HPP
#define HS_PROBLEM_002_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_002 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[1]-x[0]*x[0];
    ScalarT b = 1-x[0];
    return 100*a*a+b*b;
  }
};
}

template<class Real> 
class Problem_002 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>   NP;
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;

public:

  Problem_002() : NP( dimension_x() ) {
    NP::setLower(1,1.5);
  }

  int dimension_x() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_002::Obj>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-2.0,1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(909);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0.0504261879);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real a = std::sqrt(598)/std::sqrt(1200);
    Real b = 400*a*a*a;
    Real x[] = {2*a*std::cos(std::acos(1.0/b)/3.0),1.5};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_002_HPP
