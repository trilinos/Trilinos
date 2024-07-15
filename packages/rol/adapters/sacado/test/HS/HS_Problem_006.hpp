// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_006_HPP
#define HS_PROBLEM_006_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_006 {
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
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = 10.0*(x[1]-x[0]*x[0]);    
  }
};
}


template<class Real> 
class Problem_006 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>   NP;
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;

public:

  Problem_006() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x()  { return 2; }
  int dimension_ce() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_006::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_006::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-1.2,1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(4.84);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {1.0,1.0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_006_HPP
