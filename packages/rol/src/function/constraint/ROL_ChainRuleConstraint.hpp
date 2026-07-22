// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CHAIN_RULE_CONSTRAINT_HPP
#define ROL_CHAIN_RULE_CONSTRAINT_HPP

#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::ChainRuleConstraint
    \brief Defines a constaint formed through function composition \f$c(x)=c_o(c_i(x))\f$

    \f$c_i\mathcal{X}\to\mathcal{Y}\f$ and \f$c_o:\mathcal{Y}\to\mathcal{Z}\f$

    It is assumed that both $c_i$ and $c_o$ both of ROL::Constraint type, are
    twice differentiable and that that the range of $c_i$ is in the domain of $c_o$.
---
*/


namespace ROL {

template<typename Real>
class ChainRuleConstraint: public Constraint<Real> {
public:

  ChainRuleConstraint( const Ptr<Constraint<Real>>& outer_con,
                       const Ptr<Constraint<Real>>& inner_con,
                       const Vector<Real>&          x,
                       const Vector<Real>&          lag_inner )
  : outer_con_(outer_con), inner_con_(inner_con), 
  y_(lag_inner.dual().clone()),    // Inner constraint space vector
  Jiv_(lag_inner.dual().clone()),  // Inner Jacobian applied to optimization space vector
  aJol_(lag_inner.dual().clone()), // Outer adjoint Jacobian applied to outer dual constaint space vector
  HiaJol_(x.dual().clone()),       // Inner Hessian applied to dual optimization space vector
  HolJiv_(lag_inner.clone()),      // Outer Hessian applied to inner constaint space vector
  tol_(0) {
  }

  virtual ~ChainRuleConstraint() = default;

  virtual void update( const Vector<Real>& x,
                             UpdateType    type,
                             int           iter = -1 ) {
    inner_con_->update(x,type,iter);
    inner_con_->value(*y_,x,tol_);
    outer_con_->update(*y_,type,iter);
  }

  virtual void update( const Vector<Real>& x,
                             bool          flag,
                             int           iter = -1 ) {
    inner_con_->update(x,flag,iter);
    inner_con_->value(*y_,x,tol_);
    outer_con_->update(*y_,flag,iter);
  }

  virtual void value(       Vector<Real>& c,
                      const Vector<Real>& x, 
                            Real&         tol ) override {
    inner_con_->value(*y_,x,tol);
    outer_con_->value(c,*y_,tol);
  }

  virtual void applyJacobian(       Vector<Real>& jv,
                              const Vector<Real>& v, 
                              const Vector<Real>& x, 
                                    Real&         tol ) override {
    inner_con_->value(*y_,x,tol);
    inner_con_->applyJacobian(*Jiv_,v,x,tol);
    outer_con_->applyJacobian(jv,*Jiv_,*y_,tol);
  }

  virtual void applyAdjointJacobian(       Vector<Real>& ajl,
                                     const Vector<Real>& l, 
                                     const Vector<Real>& x, 
                                           Real&         tol ) override {
    inner_con_->value(*y_,x,tol);
    outer_con_->applyAdjointJacobian(*aJol_,l,*y_,tol);
    inner_con_->applyAdjointJacobian(ajl,*aJol_,x,tol);
  }

  virtual void applyAdjointHessian(       Vector<Real>& ahlv,
                                    const Vector<Real>& l, 
                                    const Vector<Real>& v, 
                                    const Vector<Real>& x, 
                                          Real&         tol ) override {
    inner_con_->value(*y_,x,tol);
    inner_con_->applyJacobian(*Jiv_,v,x,tol);
    outer_con_->applyAdjointJacobian(*aJol_,l,*y_,tol);
    inner_con_->applyAdjointHessian(*HiaJol_,*aJol_,v,x,tol);
    outer_con_->applyAdjointHessian(*HolJiv_,l,*Jiv_,*y_,tol);
    inner_con_->applyAdjointJacobian(ahlv,*HolJiv_,x,tol);
    ahlv.plus(*HiaJol_);
  }

private:

  const Ptr<Constraint<Real>> outer_con_, inner_con_;
  
  Ptr<Vector<Real>> y_;         // Inner constraint space vector
  Ptr<Vector<Real>> Jiv_;       // Inner Jacobian applied to optimization space vector
  Ptr<Vector<Real>> aJol_;      // Outer adjoint Jacobian applied to outer dual constaint space vector
  Ptr<Vector<Real>> HiaJol_;    // Inner Hessian applied to dual optimization space vector
  Ptr<Vector<Real>> HolJiv_;    // Outer Hessian applied to inner constaint space vector
  Real tol_;
}; // class ChainRuleConstraint

} // namespace ROL

#endif // ROL_CHAIN_RULE_CONSTRAINT_HPP
