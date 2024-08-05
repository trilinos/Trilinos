// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_020_HPP
#define HS_PROBLEM_020_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_020 {
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

template<class Real> 
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = x[0] + x[1]*x[1];
    c[1] = x[0]*x[0] + x[1];
    c[2] = x[0]*x[0] + x[1]*x[1] - 1;
  }
};
}


template<class Real> 
class Problem_020 : public ROL::NonlinearProgram<Real> {

  
  
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;
  typedef ROL::NonlinearProgram<Real>     NP;

public:

  Problem_020() : NP( dimension_x() ) {
    NP::setLower(0,-0.5);
    NP::setUpper(0,0.5);
  }  

  int dimension_x()  { return 2; }
  int dimension_ci() { return 3; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_020::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_020::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-2,1};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(909);
  }
 
  Real getSolutionObjectiveValue() {
    return 81.5 - 25*std::sqrt(3);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {0.5,0.5*std::sqrt(3)};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_020_HPP
