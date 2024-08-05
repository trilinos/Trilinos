// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_040_HPP
#define HS_PROBLEM_040_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_040 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return -x[0]*x[1]*x[2]*x[3];
  }
};

template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = x[0]*x[0]*x[0] + x[1]*x[1] - 1.0;    
    c[1] = x[0]*x[0]*x[3] - x[2];
    c[2] = x[3]*x[3]      - x[1];
  }
};
} // HS_040


template<class Real> 
class Problem_040 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>   NP;
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;

public:

  Problem_040() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x()  { return 4; }
  int dimension_ce() { return 3; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_040::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_040::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0.8,0.8,0.8,0.8};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(-0.4096);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-0.25);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real a = -1.0/3.0;
    Real b = -1.0/4.0;
    Real c = -11.0/12.0;
    Real x1[] = {std::pow(2.0,a),std::pow(2.0,2*b),-std::pow(2.0,c),-std::pow(2.0,b)};
    Real x2[] = {std::pow(2.0,a),std::pow(2.0,2*b), std::pow(2.0,c), std::pow(2.0,b)};

    return ROL::CreatePartitionedVector(NP::createOptVector(x1), 
                                        NP::createOptVector(x2));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_040_HPP
