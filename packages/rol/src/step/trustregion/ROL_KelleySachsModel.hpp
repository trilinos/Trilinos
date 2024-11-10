// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_KELLEYSACHSMODEL_HPP
#define ROL_KELLEYSACHSMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::KelleySachsModel
    \brief Provides the interface to evaluate projected trust-region model
    functions from the Kelley-Sachs bound constrained trust-region algorithm.

    -----
*/

namespace ROL {

template<class Real>
class KelleySachsModel : public TrustRegionModel<Real> {
private:
  Ptr<Vector<Real>> dual_, prim_, prim2_;
  Real eps_;

  class UpperBinding : public Elementwise::BinaryFunction<Real> {
    public:
    UpperBinding(Real offset) : offset_(offset) {}
    Real apply( const Real &x, const Real &y ) const {
      const Real one(1), tol(1e2*ROL_EPSILON<Real>());
      return ((y <= -offset_ && x <= tol) ? -one : one);
    }
    private:
    Real offset_;
  };

  class LowerBinding : public Elementwise::BinaryFunction<Real> {
    public:
    LowerBinding(Real offset) : offset_(offset) {}
    Real apply( const Real &x, const Real &y ) const {
      const Real one(1), tol(1e2*ROL_EPSILON<Real>());
      return ((y >= offset_ && x <= tol) ? -one : one);
    }
    private:
    Real offset_;
  };

  class PruneBinding : public Elementwise::BinaryFunction<Real> {
    public:
    Real apply( const Real &x, const Real &y ) const {
      const Real zero(0), one(1);
      return ((y == one) ? x : zero);
    }
  } binding_;

  class PruneNonbinding : public Elementwise::BinaryFunction<Real> {
    public:
    Real apply( const Real &x, const Real &y ) const {
      const Real zero(0), one(1);
      return ((y == -one) ? x : zero);
    }
  } nonbinding_;

  void pruneBindingConstraints(Vector<Real> &v) {
    const Ptr<const Vector<Real>> gc = TrustRegionModel<Real>::getGradient();
    const Ptr<const Vector<Real>> xc = TrustRegionModel<Real>::getIterate();
    //const Real one(1);
    //if (TrustRegionModel<Real>::getBoundConstraint()->isLowerActivated()) {
    //  prim2_->set(*xc);
    //  prim2_->axpy(-one,*TrustRegionModel<Real>::getBoundConstraint()->getLowerBound());
    //  LowerBinding op(eps_);
    //  prim2_->applyBinary(op,*gc);
    //  v.applyBinary(binding_,*prim2_);
    //}
    //if (TrustRegionModel<Real>::getBoundConstraint()->isUpperActivated()) {
    //  prim2_->set(*TrustRegionModel<Real>::getBoundConstraint()->getUpperBound());
    //  prim2_->axpy(-one,*xc);
    //  UpperBinding op(eps_);
    //  prim2_->applyBinary(op,*gc);
    //  v.applyBinary(binding_,*prim2_);
    //}
    TrustRegionModel<Real>::getBoundConstraint()->pruneActive(v,*gc,*xc,eps_);
  }

  void pruneNonbindingConstraints(Vector<Real> &v) {
    const Ptr<const Vector<Real>> gc = TrustRegionModel<Real>::getGradient();
    const Ptr<const Vector<Real>> xc = TrustRegionModel<Real>::getIterate();
    //const Real one(1);
    //if (TrustRegionModel<Real>::getBoundConstraint()->isLowerActivated()) {
    //  prim2_->set(*xc);
    //  prim2_->axpy(-one,*TrustRegionModel<Real>::getBoundConstraint()->getLowerBound());
    //  LowerBinding op(eps_);
    //  prim2_->applyBinary(op,*gc);
    //  v.applyBinary(nonbinding_,*prim2_);
    //}
    //if (TrustRegionModel<Real>::getBoundConstraint()->isUpperActivated()) {
    //  prim2_->set(*TrustRegionModel<Real>::getBoundConstraint()->getUpperBound());
    //  prim2_->axpy(-one,*xc);
    //  UpperBinding op(eps_);
    //  prim2_->applyBinary(op,*gc);
    //  v.applyBinary(nonbinding_,*prim2_);
    //}
    TrustRegionModel<Real>::getBoundConstraint()->pruneInactive(v,*gc,*xc,eps_);
  }

public:

  KelleySachsModel(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                   const Vector<Real> &x, const Vector<Real> &g,
                   const Ptr<Secant<Real>> &secant = nullPtr,
                   const bool useSecantPrecond = false, const bool useSecantHessVec = false)
    : TrustRegionModel<Real>::TrustRegionModel(obj,bnd,x,g,secant,useSecantPrecond,useSecantHessVec),
      eps_(1) {
    prim_ = x.clone();
    dual_ = g.clone();
    prim2_ = x.clone();
  }

  void setEpsilon(const Real eps) {
    eps_ = eps;
  }

  /***************************************************************************/
  /*********  BEGIN OBJECTIVE FUNCTION DEFINITIONS  **************************/
  /***************************************************************************/
  Real value( const Vector<Real> &s, Real &tol ) {
    hessVec(*dual_,s,s,tol);
    dual_->scale(static_cast<Real>(0.5));
    // Remove active components of gradient
    prim_->set(TrustRegionModel<Real>::getGradient()->dual());
    pruneBindingConstraints(*prim_);
    // Add reduced gradient to reduced hessian in direction s
    dual_->plus(prim_->dual());
    return dual_->dot(s.dual());
  }

  void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) {
    // Apply (reduced) hessian to direction s
    hessVec(g,s,s,tol);
    // Remove active components of gradient
    prim_->set(TrustRegionModel<Real>::getGradient()->dual());
    pruneBindingConstraints(*prim_);
    // Add reduced gradient to reduced hessian in direction s
    g.plus(prim_->dual());
  }

  void hessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    // Set vnew to v
    prim_->set(v);
    // Remove elements of vnew corresponding to binding set
    pruneBindingConstraints(*prim_);
    // Apply full Hessian to reduced vector
    TrustRegionModel<Real>::applyHessian(Hv,*prim_,tol);
    // Remove elements of Hv corresponding to binding set
    pruneBindingConstraints(Hv);
    // Set vnew to v
    prim_->set(v);
    // Remove Elements of vnew corresponding to complement of binding set
    pruneNonbindingConstraints(*prim_);
    dual_->set(prim_->dual());
    pruneNonbindingConstraints(*dual_);
    // Fill complement of binding set elements in Hp with v
    Hv.plus(*dual_);
  }

  void invHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    // Set vnew to v
    dual_->set(v);
    // Remove elements of vnew corresponding to binding set
    pruneBindingConstraints(*dual_);
    // Apply full Hessian to reduced vector
    TrustRegionModel<Real>::applyInvHessian(Hv,*dual_,tol);
    // Remove elements of Hv corresponding to binding set
    pruneBindingConstraints(Hv);
    // Set vnew to v
    dual_->set(v);
    // Remove Elements of vnew corresponding to complement of binding set
    pruneNonbindingConstraints(*dual_);
    prim_->set(dual_->dual());
    pruneNonbindingConstraints(*prim_);
    // Fill complement of binding set elements in Hv with v
    Hv.plus(*prim_);
  }

  void precond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Set vnew to v
    dual_->set(v);
    // Remove elements of vnew corresponding to binding set
    pruneBindingConstraints(*dual_);
    // Apply full Hessian to reduced vector
    TrustRegionModel<Real>::applyPrecond(Mv,*dual_,tol);
    // Remove elements of Mv corresponding to binding set
    pruneBindingConstraints(Mv);
    // Set vnew to v
    dual_->set(v);
    // Remove Elements of vnew corresponding to complement of binding set
    pruneNonbindingConstraints(*dual_);
    prim_->set(dual_->dual());
    pruneNonbindingConstraints(*prim_);
    // Fill complement of binding set elements in Mv with v
    Mv.plus(*prim_);
  }
  /***************************************************************************/
  /*********  END OBJECTIVE FUNCTION DEFINITIONS  ****************************/
  /***************************************************************************/

  void dualTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    // Compute T(v) = P_I(v) where P_I is the projection onto the inactive indices
    tv.set(v);
    pruneBindingConstraints(tv);
  }

  void primalTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    // Compute T(v) = P( x + v ) - x where P is the projection onto the feasible set
    const Ptr<const Vector<Real>> xc = TrustRegionModel<Real>::getIterate();
    tv.set(*xc);
    tv.plus(v);
    TrustRegionModel<Real>::getBoundConstraint()->project(tv);
    tv.axpy(static_cast<Real>(-1),*xc);
  }

}; 

} // namespace ROL

#endif 
