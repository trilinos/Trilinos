// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_019_HPP
#define HS_PROBLEM_019_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_019 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
     ScalarT a = x[0]-10;
     ScalarT b = x[1]-20;
 
    return a*a*a + b*b*b;
  }
};

template<class Real> 
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    ScalarT a = x[0]-5;
    ScalarT b = x[1]-5;
    ScalarT d = x[0]-6;

    c[0] =  a*a + b*b - 100;
    c[1] = -b*b - d*d + 82.81;
  }
};
}


template<class Real> 
class Problem_019 : public ROL::NonlinearProgram<Real> {

  
  
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;
  typedef ROL::NonlinearProgram<Real>     NP;

public:

  Problem_019() : NP( dimension_x() ) {
    NP::setLower(0,13.0);
    NP::setUpper(0,100.0);
    NP::setLower(1,0.0);
    NP::setUpper(1,100.0);
  }  

  int dimension_x()  { return 2; }
  int dimension_ci() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_019::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_019::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {20.1,5.84};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(-1808.858296);
  }
 
  Real getSolutionObjectiveValue() {
    return -6961.81381;
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {14.095,0.84296079};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_019_HPP
