// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_INTERIORPOINT_PRIMALDUAL_RESIDUAL_H
#define ROL_INTERIORPOINT_PRIMALDUAL_RESIDUAL_H

#include "ROL_Elementwise_Function.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_Objective.hpp"
#include "ROL_PartitionedVector.hpp"

namespace ROL {
namespace InteriorPoint {

/** @ingroup func_group
    \class   ROL::InteriorPoint::PrimalDualResidual
    \brief   Express the Primal-Dual Interior Point gradient
             as an equality constraint 
    
             See Nocedal & Wright second edition equation (19.6)
             In that book the convention for naming components
      
             x - optimization variable (here subscript o)
             s - slack variable (here subscript s)
             y - Lagrange multiplier for the equality constraint (here subscript e)
             z - Lagrange multiplier for the inequality constraint (here subscript i)
    --- 
*/


template<class Real> class PrimalDualSymmetrizer;


template<class Real>
class PrimalDualResidual : public Constraint<Real> {

private:
  typedef Vector<Real>             V;
  typedef PartitionedVector<Real>  PV;
  typedef Objective<Real>          OBJ;
  typedef Constraint<Real> CON;


  typedef typename PV::size_type   size_type;

  ROL::Ptr<OBJ> obj_;    // Objective function
  ROL::Ptr<CON> eqcon_;  //  Constraint
  ROL::Ptr<CON> incon_;  // Inequality Constraint

  ROL::Ptr<V> qo_;       // Storage for optimization variables
  ROL::Ptr<V> qs_;       // Storage for slack variables
  ROL::Ptr<V> qe_;       // Storage for equality multiplier variables
  ROL::Ptr<V> qi_;       // Storage for inequality multiplier variables

  Real mu_;                  // Penalty parameter

  ROL::Ptr<LinearOperator<Real> > sym_;

  const static size_type OPT   = 0;  // Optimization vector
  const static size_type SLACK = 1;  // Slack vector
  const static size_type EQUAL = 2;  // Lagrange multipliers for equality constraint
  const static size_type INEQ  = 3;  // Lagrange multipliers for inequality constraint



public:

  PrimalDualResidual( const ROL::Ptr<OBJ> &obj, 
                      const ROL::Ptr<CON> &eqcon,
                      const ROL::Ptr<CON> &incon,
                      const V& x ) : 
                      obj_(obj), eqcon_(eqcon), incon_(incon), mu_(1.0) {

    // Allocate storage vectors
    const PV &xpv = dynamic_cast<const PV&>(x);

    qo_ = xpv.get(OPT)->clone();
    qs_ = xpv.get(SLACK)->clone();
    qe_ = xpv.get(EQUAL)->clone();
    qi_ = xpv.get(INEQ)->clone();

    sym_ = ROL::makePtr<PrimalDualSymmetrizer<Real>>(*qs_);

  }

  void value( V &c, const V &x, Real &tol ) {

    

    // Downcast to partitioned vectors
    PV &cpv = dynamic_cast<PV&>(c);
    const PV &xpv = dynamic_cast<const PV&>(x);

    ROL::Ptr<const V> xo = xpv.get(OPT);
    ROL::Ptr<const V> xs = xpv.get(SLACK);
    ROL::Ptr<const V> xe = xpv.get(EQUAL);
    ROL::Ptr<const V> xi = xpv.get(INEQ);

    c.zero();    

    ROL::Ptr<V> co = cpv.get(OPT);
    ROL::Ptr<V> cs = cpv.get(SLACK);
    ROL::Ptr<V> ce = cpv.get(EQUAL);
    ROL::Ptr<V> ci = cpv.get(INEQ);

    // Optimization components
    obj_->gradient(*co,*xo,tol); 

    // Apply adjoint equality Jacobian at xo to xe and store in qo
    eqcon_->applyAdjointJacobian(*qo_,*xe,*xo,tol);
    co->axpy(-1.0,*qo_);

    incon_->applyAdjointJacobian(*qo_,*xi,*xo,tol);
    co->axpy(-1.0,*qo_);

    // Slack components
    cs->set(*xs);

    Elementwise::Multiply<Real> mult;
    cs->applyBinary(mult,*xi);
 
    Elementwise::Fill<Real> fill(-mu_);
    qs_->applyUnary(fill);

    cs->plus(*qs_);               // Sz-e

    //  component
    eqcon_->value(*ce, *xo, tol); // c_e(x)

    // Inequality component
    incon_->value(*ci, *xo, tol); // c_i(x)-s
    ci->axpy(-1.0, *xs); 

    sym_->update(*xs);
    sym_->apply(c,c,tol);
    sym_->applyInverse(c,c,tol);

  }   

  void applyJacobian( V &jv, const V &v, const V &x, Real &tol ) {

    

    jv.zero();

    // Downcast to partitioned vectors
    PV &jvpv = dynamic_cast<PV&>(jv);
    const PV &vpv = dynamic_cast<const PV&>(v);
    const PV &xpv = dynamic_cast<const PV&>(x);

    ROL::Ptr<V> jvo = jvpv.get(OPT);
    ROL::Ptr<V> jvs = jvpv.get(SLACK);
    ROL::Ptr<V> jve = jvpv.get(EQUAL);
    ROL::Ptr<V> jvi = jvpv.get(INEQ);

    ROL::Ptr<const V> vo = vpv.get(OPT);
    ROL::Ptr<const V> vs = vpv.get(SLACK);
    ROL::Ptr<const V> ve = vpv.get(EQUAL);
    ROL::Ptr<const V> vi = vpv.get(INEQ);

    ROL::Ptr<const V> xo = xpv.get(OPT);
    ROL::Ptr<const V> xs = xpv.get(SLACK);
    ROL::Ptr<const V> xe = xpv.get(EQUAL);
    ROL::Ptr<const V> xi = xpv.get(INEQ);

    // Optimization components
    obj_->hessVec(*jvo,*vo,*xo,tol);
    
    eqcon_->applyAdjointHessian(*qo_,*xe,*vo,*xo,tol);

    jvo->axpy(-1.0,*qo_);

    incon_->applyAdjointHessian(*qo_,*xi,*vo,*xo,tol);
    
    jvo->axpy(-1.0,*qo_);
    
    eqcon_->applyAdjointJacobian(*qo_,*ve,*xo,tol);

    jvo->axpy(-1.0,*qo_);
 
    incon_->applyAdjointJacobian(*qo_,*vi,*xo,tol);

    jvo->axpy(-1.0,*qo_);
    

    // Slack components
    jvs->set(*vs);

    Elementwise::Multiply<Real> mult;
 
    jvs->applyBinary(mult,*xi);

    qs_->set(*vi);

    qs_->applyBinary(mult,*xs);
    
    jvs->plus(*qs_);

    //  component
    eqcon_->applyJacobian(*jve,*vo,*xo,tol);
    
    // Inequality components
    incon_->applyJacobian(*jvi,*vo,*xo,tol);
    
    jvi->axpy(-1.0,*vs);

    sym_->update(*xs);
    sym_->apply(jv,jv,tol);
    sym_->applyInverse(jv,jv,tol); 



  }

  void updatePenalty( Real mu ) { 
    mu_ = mu;
  }

};  // class PrimalDualResidual



// Applying this operator to the left- and right-hand-sides, will
// symmetrize the Primal-Dual Interior-Point KKT system, yielding
// equation (19.13) from Nocedal & Wright

template<class Real> 
class PrimalDualSymmetrizer : public LinearOperator<Real> {

  typedef Vector<Real>             V;  
  typedef PartitionedVector<Real>  PV;

  typedef typename PV::size_type   size_type;

private:
  ROL::Ptr<V> s_;

  const static size_type OPT   = 0;  // Optimization vector
  const static size_type SLACK = 1;  // Slack vector
  const static size_type EQUAL = 2;  // Lagrange multipliers for equality constraint
  const static size_type INEQ  = 3;  // Lagrange multipliers for inequality constraint

public:
  
  PrimalDualSymmetrizer( const V &s ) {
    s_ = s.clone();
    s_->set(s);
  }

  void update( const V& s, bool flag = true, int iter = -1 ) {
    s_->set(s);  
  }

  void apply( V &Hv, const V &v, Real &tol ) const {

    
    

    const PV &vpv = dynamic_cast<const PV&>(v);
    PV &Hvpv = dynamic_cast<PV&>(Hv);

    ROL::Ptr<const V> vo = vpv.get(OPT);
    ROL::Ptr<const V> vs = vpv.get(SLACK);
    ROL::Ptr<const V> ve = vpv.get(EQUAL);
    ROL::Ptr<const V> vi = vpv.get(INEQ);

    ROL::Ptr<V> Hvo = Hvpv.get(OPT);
    ROL::Ptr<V> Hvs = Hvpv.get(SLACK);
    ROL::Ptr<V> Hve = Hvpv.get(EQUAL);
    ROL::Ptr<V> Hvi = Hvpv.get(INEQ);

    Hvo->set(*vo);

    Hvs->set(*vs);
    Elementwise::Divide<Real> div;
    Hvs->applyBinary(div,*s_);

    Hve->set(*ve);
    Hve->scale(-1.0);

    Hvi->set(*vi);
    Hvi->scale(-1.0);

  }

  void applyInverse( V &Hv, const V &v, Real &tol ) const {

    
    

    const PV &vpv = dynamic_cast<const PV&>(v);
    PV &Hvpv = dynamic_cast<PV&>(Hv);

    ROL::Ptr<const V> vo = vpv.get(OPT);
    ROL::Ptr<const V> vs = vpv.get(SLACK);
    ROL::Ptr<const V> ve = vpv.get(EQUAL);
    ROL::Ptr<const V> vi = vpv.get(INEQ);

    ROL::Ptr<V> Hvo = Hvpv.get(OPT);
    ROL::Ptr<V> Hvs = Hvpv.get(SLACK);
    ROL::Ptr<V> Hve = Hvpv.get(EQUAL);
    ROL::Ptr<V> Hvi = Hvpv.get(INEQ);

    Hvo->set(*vo);

    Hvs->set(*vs);
    Elementwise::Multiply<Real> mult;
    Hvs->applyBinary(mult,*s_);

    Hve->set(*ve);
    Hve->scale(-1.0);

    Hvi->set(*vi);
    Hvi->scale(-1.0);

  }
}; // class PrimalDualSymmetrizer


} // namespace InteriorPoint





} // namespace ROL


#endif // ROL_INTERIORPOINT_PRIMALDUAL_RESIDUAL_H

