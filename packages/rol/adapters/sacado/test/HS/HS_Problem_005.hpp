// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_005_HPP
#define HS_PROBLEM_005_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_005 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]+x[1];
    ScalarT b = x[1]-x[0];
    return std::sin(a) + b*b - 1.5*x[0] + 2.5*x[1] + 1.0;
  }
};
}


template<class Real> 
class Problem_005 : public ROL::NonlinearProgram<Real> {


  

  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::NonlinearProgram<Real>   NP;

  const Real pi;

public:

  Problem_005() : NP( dimension_x() ), pi(3.14159265459) {
    Real l[] = {-1.5,-3.0};
    Real u[] = {4.0,3.0};
    NP::setLower(l); NP::setUpper(u);
  }

  int dimension_x() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_005::Obj>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0,0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(1.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-std::sqrt(3.0)/2.0-pi/3.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {-pi/3.0+0.5,-pi/3.0-0.5};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_005_HPP
