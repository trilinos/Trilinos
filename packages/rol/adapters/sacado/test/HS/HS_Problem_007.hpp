// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_007_HPP
#define HS_PROBLEM_007_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_007 {

template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = 1.0+x[0]*x[0];
    return std::log(a) - x[1];
  }
};

template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    ScalarT a = 1.0 + x[0]*x[0];
    c[0] = a*a + x[1]*x[1] - 4;    
  }
};

}

template<class Real> 
class Problem_007 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;
  typedef ROL::NonlinearProgram<Real>   NP;

public:

  Problem_007() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x()  { return 2; }
  int dimension_ce() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_007::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_007::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {2,2};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return std::log(5.0)-2.0;
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-std::sqrt(3.0));
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {0.0, std::sqrt(3)};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_007_HPP
