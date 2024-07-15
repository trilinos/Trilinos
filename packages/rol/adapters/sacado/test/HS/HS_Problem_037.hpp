// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_037_HPP
#define HS_PROBLEM_037_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_037 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return -x[0]*x[1]*x[2];
 }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {

    c[0] = 72.0 - x[0] - 2.0*x[1] - 2.0*x[2];
    c[1] =        x[0] + 2.0*x[1] + 2.0*x[2];

  }
};
} // namespace HS_037


template<class Real> 
class Problem_037 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

public:

  Problem_037() : NP( dimension_x() ) {
    NP::setLower(0,  0.0);
    NP::setUpper(0, 42.0);
    NP::setLower(1,  0.0);
    NP::setUpper(1, 42.0);
    NP::setLower(2,  0.0);  
    NP::setUpper(2, 42.0);
  }

  int dimension_x()  { return 3; }
  int dimension_ci() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_037::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_037::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {10.0,10.0,10.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(-1000.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-3456.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = { 24.0, 12.0, 12.0 };

    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_037_HPP
