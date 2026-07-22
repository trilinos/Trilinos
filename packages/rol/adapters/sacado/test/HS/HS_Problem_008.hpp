// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_008_HPP
#define HS_PROBLEM_008_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_008 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return ScalarT(-1.0);
  }
};

template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
  c[0] = x[0]*x[0] + x[1]*x[1] - 25.0;    
  c[1] = x[0]*x[1] - 9.0;
  }
};
}


template<class Real> 
class Problem_008 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::Vector<Real>              V;
  typedef ROL::Objective<Real>           OBJ;
  typedef ROL::Constraint<Real>          CON;
  typedef ROL::NonlinearProgram<Real>    NP;

public:

  Problem_008() : NP( dimension_x() ) {
    NP::noBound(); 
  }  

  int dimension_x()  { return 2; }
  int dimension_ce() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_008::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_008::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {2,1};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(-1.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-1.0);
  }

  ROL::Ptr<const V> getSolutionSet() {

    Real a = std::sqrt( (25.0 + std::sqrt(301.0))/2.0 );
    Real b = std::sqrt( (25.0 - std::sqrt(301.0))/2.0 );

    Real x1[] = { a,  9/a};
    Real x2[] = {-a, -9/a};
    Real x3[] = { b,  9/b};
    Real x4[] = {-b, -9/b};

    return ROL::CreatePartitionedVector(NP::createOptVector(x1),
                                        NP::createOptVector(x2),
                                        NP::createOptVector(x3),
                                        NP::createOptVector(x4));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_008_HPP
