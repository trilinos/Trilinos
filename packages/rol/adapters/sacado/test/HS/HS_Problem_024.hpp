// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_024_HPP
#define HS_PROBLEM_024_HPP

#include "ROL_NonlinearProgram.hpp"


namespace HS {

namespace HS_024 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]-3;
    return (a*a-9)*x[1]*x[1]*x[1]/(27*std::sqrt(3));
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
    void value( std::vector<ScalarT> &c,
                const std::vector<ScalarT> &x,
                Real &tol ) {
    Real q = std::sqrt(3);
    c[0] =  x[0]/q -   x[1];
    c[1] =  x[0]   + q*x[1];
    c[2] = -x[0]   - q*x[1] + 6;
  }
};
}




template<class Real> 
class Problem_024 : public ROL::NonlinearProgram<Real> {
 
  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

private:
public:

  Problem_024() : NP( dimension_x() ) {
    NP::setLower(0,0.0);
    NP::setLower(1,0.0);
  }

  int dimension_x() { return 2; }
  int dimension_ci() { return 3; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_024::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_024::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {1.0,0.5};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(-0.01336459);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-1.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    const Real x[] = {3,std::sqrt(3)};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

}

#endif // HS_PROBLEM_024_HPP
