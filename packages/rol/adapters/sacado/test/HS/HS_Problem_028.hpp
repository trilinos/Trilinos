// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_028_HPP
#define HS_PROBLEM_028_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_028 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = (x[0]+x[1]);
    ScalarT b = (x[1]+x[2]);
    return a*a + b*b;
  }
};


template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = x[0]+2*x[1]+3*x[2]-1.0;
  }
};
}


template<class Real> 
class Problem_028 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>   NP;
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;

public:

  Problem_028() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x()  { return 3; }
  int dimension_ce() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_028::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_028::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-4.0,1.0,1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(13.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {0.5,-0.5,0.5};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_028_HPP
