// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_035_HPP
#define HS_PROBLEM_035_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_035 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return 9.0 - 8.0*x[0]      - 6.0*x[1]      - 4.0*x[2]
               + 2.0*x[0]*x[0] + 2.0*x[1]*x[1] +     x[2]*x[2] 
               + 2.0*x[0]*x[1] + 2.0*x[0]*x[2];
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {

    c[0] = 3.0 - x[0] - x[1] - 2.0*x[2];

  }
};
} // namspace HS_035


template<class Real> 
class Problem_035 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

public:

  Problem_035() : NP( dimension_x() ) {
    NP::setLower(0, 0.0);
    NP::setLower(1, 0.0);
    NP::setLower(2, 0.0);
  }

  int dimension_x()  { return 3; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_035::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_035::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0.5,0.5,0.5};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(2.25);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(1.0/9.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = { 4.0/3.0, 7.0/9.0, 4.0/9.0 };

    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_035_HPP
