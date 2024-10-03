// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MOREAUYOSIDAOBJECTIVE_H
#define ROL_MOREAUYOSIDAOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_ScalarController.hpp"
#include "ROL_ParameterList.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::MoreauYosidaObjective
    \brief Provides the interface to evaluate the Moreau-Yosida penalty function.

    ---
*/


namespace ROL {

template <class Real>
class MoreauYosidaObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  const Ptr<BoundConstraint<Real>> bnd_;

  Ptr<Vector<Real>> l_;
  Ptr<Vector<Real>> u_;
  Ptr<Vector<Real>> l1_;
  Ptr<Vector<Real>> u1_;
  Ptr<Vector<Real>> dl1_;
  Ptr<Vector<Real>> du1_;
  Ptr<Vector<Real>> xlam_;
  Ptr<Vector<Real>> v_;
  Ptr<Vector<Real>> dv_;
  Ptr<Vector<Real>> dv2_;
  Ptr<Vector<Real>> lam_;
  Ptr<Vector<Real>> tmp_;

  Ptr<ScalarController<Real,int>> fval_;
  Ptr<VectorController<Real,int>> gradient_;

  Real mu_;
  bool isPenEvaluated_;
  int nfval_;
  int ngrad_;
  bool updateMultiplier_;
  bool updatePenalty_;

  void computePenalty(const Vector<Real> &x) {
    if ( bnd_->isActivated() ) {
      Real one = 1.0;
      if ( !isPenEvaluated_ ) {
        xlam_->set(x);
        xlam_->axpy(one/mu_,*lam_);

        if ( bnd_->isFeasible(*xlam_) ) {
          l1_->zero(); dl1_->zero();
          u1_->zero(); du1_->zero();
        }
        else {
          // Compute lower penalty component
          l1_->set(*l_);
          bnd_->pruneLowerInactive(*l1_,*xlam_);
          tmp_->set(*xlam_);
          bnd_->pruneLowerInactive(*tmp_,*xlam_);
          l1_->axpy(-one,*tmp_);

          // Compute upper penalty component
          u1_->set(*xlam_);
          bnd_->pruneUpperInactive(*u1_,*xlam_);
          tmp_->set(*u_);
          bnd_->pruneUpperInactive(*tmp_,*xlam_);
          u1_->axpy(-one,*tmp_);

          // Compute derivative of lower penalty component
          dl1_->set(l1_->dual());
          bnd_->pruneLowerInactive(*dl1_,*xlam_);

          // Compute derivative of upper penalty component
          du1_->set(u1_->dual());
          bnd_->pruneUpperInactive(*du1_,*xlam_);
        }

        isPenEvaluated_ = true;
      }
    }
  }

  void initialize(const Vector<Real> &x,
                  const Vector<Real> &g) {
    fval_     = makePtr<ScalarController<Real,int>>();
    gradient_ = makePtr<VectorController<Real,int>>();

    l_    = x.clone();
    l1_   = x.clone();
    dl1_  = g.clone();
    u_    = x.clone();
    u1_   = x.clone();
    du1_  = g.clone();
    xlam_ = x.clone();
    v_    = x.clone();
    dv_   = g.clone();
    dv2_  = g.clone();
    lam_  = x.clone();
    tmp_  = x.clone();

    l_->set(*bnd_->getLowerBound());
    u_->set(*bnd_->getUpperBound());

    lam_->zero();
    //lam_->set(*u_);
    //lam_->plus(*l_);
    //lam_->scale(0.5);
  }

public:
  MoreauYosidaObjective(const Ptr<Objective<Real>> &obj,
                      const Ptr<BoundConstraint<Real>> &bnd, 
                      const Vector<Real> &x,
                      const Vector<Real> &g,
                      const Real mu = 1e1,
                      const bool updateMultiplier = true,
                      const bool updatePenalty = true)
    : obj_(obj), bnd_(bnd), mu_(mu),
      isPenEvaluated_(false), nfval_(0), ngrad_(0),
      updateMultiplier_(updateMultiplier), updatePenalty_(updatePenalty) {
    initialize(x,g);
  }

  MoreauYosidaObjective(const Ptr<Objective<Real>> &obj,
                      const Ptr<BoundConstraint<Real>> &bnd, 
                      const Vector<Real> &x,
                      const Vector<Real> &g,
                      const Vector<Real> &lam,
                      const Real mu = 1e1,
                      const bool updateMultiplier = true,
                      const bool updatePenalty = true)
    : MoreauYosidaObjective(obj,bnd,x,g,mu,updateMultiplier,updatePenalty) {
    lam_->set(lam);
  }

  MoreauYosidaObjective(const Ptr<Objective<Real>> &obj,
                      const Ptr<BoundConstraint<Real>> &bnd,
                      const Vector<Real> &x,
                      const Vector<Real> &g,
                      ParameterList &parlist)
    : obj_(obj), bnd_(bnd),
      isPenEvaluated_(false), nfval_(0), ngrad_(0) {
    initialize(x,g);
    ParameterList &list = parlist.sublist("Step").sublist("Moreau-Yosida Penalty");
    updateMultiplier_ = list.get("Update Multiplier",true);
    updatePenalty_    = list.get("Update Penalty",true);
    mu_               = list.get("Initial Penalty Parameter",1e1);
  }

  MoreauYosidaObjective(const Ptr<Objective<Real>> &obj,
                      const Ptr<BoundConstraint<Real>> &bnd,
                      const Vector<Real> &x,
                      const Vector<Real> &g,
                      const Vector<Real> &lam,
                      ParameterList &parlist)
    : MoreauYosidaObjective(obj,bnd,x,g,parlist) {
    lam_->set(lam);
  }

  void updateMultipliers(Real mu, const Vector<Real> &x) {
    if ( bnd_->isActivated() ) {
      if ( updateMultiplier_ ) {
        const Real one(1);
        computePenalty(x);
        lam_->set(*u1_);
        lam_->axpy(-one,*l1_);
        lam_->scale(mu_);
      }
      if ( updatePenalty_ ) {
        mu_ = mu;
      }
    }
    nfval_ = 0; ngrad_ = 0;
    isPenEvaluated_ = false;
  }

  void reset(const Real mu) {
    lam_->zero();
    mu_ = mu;
    nfval_ = 0; ngrad_ = 0;
    isPenEvaluated_ = false;
  }

  Real testComplementarity(const Vector<Real> &x) {
    Real val(0);
    if (bnd_->isActivated()) {
      computePenalty(x);

      tmp_->set(x);
      tmp_->axpy(static_cast<Real>(-1), *l_);
      Real lower = mu_*std::abs(tmp_->dot(*l1_));

      tmp_->set(x);
      tmp_->axpy(static_cast<Real>(-1), *u_);
      Real upper = mu_*std::abs(tmp_->dot(*u1_));

      tmp_->set(x);
      bnd_->project(*tmp_);
      tmp_->axpy(static_cast<Real>(-1), x);
      Real xnorm = tmp_->norm();

      val = std::max(xnorm,std::max(lower,upper));
    }
    return val;
  }

  Real getObjectiveValue(const Vector<Real> &x, Real &tol) {
    int key(0);
    Real val(0);
    bool isComputed = fval_->get(val,key);
    if (!isComputed) {
      val = obj_->value(x,tol); nfval_++;
      fval_->set(val,key);
    }
    return val;
  }

  void getObjectiveGradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {
    int key(0);
    bool isComputed = gradient_->get(g,key);
    if (!isComputed) {
      obj_->gradient(g,x,tol); ngrad_++;
      gradient_->set(g,key);
    }
  }

  int getNumberFunctionEvaluations(void) {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) {
    return ngrad_;
  }

  /** \brief Update Moreau-Yosida penalty function. 

      This function updates the Moreau-Yosida penalty function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    obj_->update(x,type,iter);
    fval_->objectiveUpdate(type);
    gradient_->objectiveUpdate(type);
    // Need to do something smart here
    isPenEvaluated_  = false;
  }

  /** \brief Compute value.

      This function returns the Moreau-Yosida penalty value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact Moreau-Yosida penalty computation.
  */
  Real value( const Vector<Real> &x, Real &tol ) {
    // Compute objective function value
    Real fval = getObjectiveValue(x,tol);
    // Add value of the Moreau-Yosida penalty
    if ( bnd_->isActivated() ) {
      const Real half(0.5);
      computePenalty(x);
      fval += half*mu_*(l1_->dot(*l1_) + u1_->dot(*u1_));
    }
    return fval;
  }

  /** \brief Compute gradient.

      This function returns the Moreau-Yosida penalty gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact Moreau-Yosida penalty computation.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Compute gradient of objective function
    getObjectiveGradient(g,x,tol);
    // Add gradient of the Moreau-Yosida penalty
    if ( bnd_->isActivated() ) {
      computePenalty(x);
      g.axpy(-mu_,*dl1_);
      g.axpy(mu_,*du1_);
    }
  }

  /** \brief Apply Hessian approximation to vector.

      This function applies the Hessian of the Moreau-Yosida penalty to the vector \f$v\f$.
      @param[out]         hv  is the the action of the Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact Moreau-Yosida penalty computation.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    // Add Hessian of the Moreau-Yosida penalty
    if ( bnd_->isActivated() ) {
      const Real one(1);
      computePenalty(x);

      v_->set(v);
      bnd_->pruneLowerActive(*v_,*xlam_);
      v_->scale(-one);
      v_->plus(v);
      dv_->set(v_->dual());
      dv2_->set(*dv_);
      bnd_->pruneLowerActive(*dv_,*xlam_);
      dv_->scale(-one);
      dv_->plus(*dv2_);
      hv.axpy(mu_,*dv_);

      v_->set(v);
      bnd_->pruneUpperActive(*v_,*xlam_);
      v_->scale(-one);
      v_->plus(v);
      dv_->set(v_->dual());
      dv2_->set(*dv_);
      bnd_->pruneUpperActive(*dv_,*xlam_);
      dv_->scale(-one);
      dv_->plus(*dv2_);
      hv.axpy(mu_,*dv_);
    }
  }
}; // class MoreauYosidaObjective

} // namespace ROL

#endif
