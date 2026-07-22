// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LOG_BARRIER_OBJECTIVE_H
#define ROL_LOG_BARRIER_OBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::LogBarrierObjective
    \brief Log barrier objective for interior point methods
*/

namespace ROL {

template <class Real>
class LogBarrierObjective : public Objective<Real> {
public:

    

  /* \brief Objective value J(x) = \f$-\sum_i \log(x_i) \f$ */
  Real value( const Vector<Real> &x, Real &tol ) {

    ROL::Ptr<Vector<Real> > logx = x.clone();
    logx->set(x);
    
    Elementwise::Logarithm<Real> log;

    logx->applyUnary(log);

    Elementwise::ReductionSum<Real> sum;

    Real result = -(logx->reduce(sum));

    return result;
  }  

  /* \brief gradient g_i = \f$ 1/x_i \f$ */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    
    g.set(x);
     
    Elementwise::Reciprocal<Real> reciprocal;
   
    g.applyUnary(reciprocal);
    g.scale(-1.0);
  }

  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
 
    ROL::Ptr<Vector<Real> > dbyx = d.clone();
    dbyx->set(x);

    struct Division : public Elementwise::BinaryFunction<Real> {
      Real apply( const Real &xc, const Real &dc ) const {
        return dc/xc;    
      }
    } division;

    dbyx->applyBinary( division, d );
    
    Elementwise::ReductionSum<Real> sum;

    return -dbyx->reduce(sum);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.set(v);
    
    struct HessianApply : public Elementwise::BinaryFunction<Real> {
      Real apply( const Real &vc, const Real &xc ) const { 
        return vc/(xc*xc);   
      }
    } hessian;

    hv.applyBinary(hessian,x);

  }

};

} // namespace ROL

#endif // ROL_LOG_BARRIER_OBJECTIVE_H
