// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_022_HPP
#define HS_PROBLEM_022_HPP

#include "ROL_NonlinearProgram.hpp"


namespace HS {

namespace HS_022 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]-2;
    ScalarT b = x[1]-1;
    return a*a+b*b;
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
    void value( std::vector<ScalarT> &c,
                const std::vector<ScalarT> &x,
                Real &tol ) {
    c[0] = -x[0]-x[1]+2;
    c[1] = -x[0]*x[0]+x[1];
  }
};
}




template<class Real> 
class Problem_022 : public ROL::NonlinearProgram<Real> {
 
  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

private:
public:

  Problem_022() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x() { return 2; }
  int dimension_ci() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_022::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_022::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {2.0,2.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(1.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(1);
  }

  ROL::Ptr<const V> getSolutionSet() {
    const Real x[] = {1.0,1.0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

}

#endif // HS_PROBLEM_022_HPP
