// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTEDLAGRANGIANOBJECTIVE_H
#define ROL_AUGMENTEDLAGRANGIANOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_QuadraticPenalty.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_ScalarController.hpp"
#include "ROL_ParameterList.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::AugmentedLagrangianObjective
    \brief Provides the interface to evaluate the augmented Lagrangian.

    This class implements the augmented Lagrangian functional for use with
    ROL::AugmentedLagrangianAlgorithm.  Given a function
    \f$f:\mathcal{X}\to\mathbb{R}\f$ and an equality constraint
    \f$c:\mathcal{X}\to\mathcal{C}\f$, the augmented Lagrangian functional is
    \f[
       L_A(x,\lambda,\mu) = f(x) +
           \langle \lambda, c(x)\rangle_{\mathcal{C}^*,\mathcal{C}} +
           \frac{\mu}{2} \langle \mathfrak{R}c(x),c(x)\rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$\lambda\in\mathcal{C}^*\f$ denotes the Lagrange multiplier estimate,
    \f$\mu > 0\f$ is the penalty parameter and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ is the Riesz operator
    on the constraint space.

    This implementation permits the scaling of \f$L_A\f$ by \f$\mu^{-1}\f$ and also
    permits the Hessian approximation
    \f[
        \nabla^2_x L_A(x,\lambda,\mu)v \approx \nabla^2 f(x) v + \mu c'(x)^*\mathfrak{R} c'(x)v.
    \f]

    ---
*/


namespace ROL {

template <class Real>
class AugmentedLagrangianObjective : public Objective<Real> {
private:
  // Required for Augmented Lagrangian definition
  const Ptr<Objective<Real>>  obj_;
  const Ptr<Constraint<Real>> con_;

  Real penaltyParameter_;
  Ptr<Vector<Real>> multiplier_;

  // Auxiliary storage
  Ptr<Vector<Real>> dualOptVector_;
  Ptr<Vector<Real>> dualConVector_;
  Ptr<Vector<Real>> primConVector_;

  // Objective and constraint evaluations
  Ptr<ScalarController<Real,int>> fval_;
  Ptr<VectorController<Real,int>> gradient_;
  Ptr<VectorController<Real,int>> conValue_;

  // Objective function and constraint scaling
  Real fscale_;
  Real cscale_;

  // Evaluation counters
  int nfval_;
  int ngval_;
  int ncval_;

  // User defined options
  bool scaleLagrangian_;
  int HessianApprox_;

public:
  AugmentedLagrangianObjective(const Ptr<Objective<Real>> &obj,
                               const Ptr<Constraint<Real>> &con,
                               const Real penaltyParameter,
                               const Vector<Real> &dualOptVec,
                               const Vector<Real> &primConVec,
                               const Vector<Real> &dualConVec,
                               ParameterList &parlist)
    : obj_(obj), con_(con), penaltyParameter_(penaltyParameter),
      fscale_(1), cscale_(1), nfval_(0), ngval_(0), ncval_(0) {

    fval_     = makePtr<ScalarController<Real,int>>();
    gradient_ = makePtr<VectorController<Real,int>>();
    conValue_ = makePtr<VectorController<Real,int>>();

    multiplier_    = dualConVec.clone();
    dualOptVector_ = dualOptVec.clone();
    dualConVector_ = dualConVec.clone();
    primConVector_ = primConVec.clone();

    ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_ = sublist.get("Use Scaled Augmented Lagrangian", false);
    HessianApprox_   = sublist.get("Level of Hessian Approximation",  0);
  }

  AugmentedLagrangianObjective(const Ptr<Objective<Real>> &obj,
                               const Ptr<Constraint<Real>> &con,
                               const Real penaltyParameter,
                               const Vector<Real> &dualOptVec,
                               const Vector<Real> &primConVec,
                               const Vector<Real> &dualConVec,
                               const bool scaleLagrangian,
                               const int HessianApprox)
    : obj_(obj), con_(con), penaltyParameter_(penaltyParameter),
      fscale_(1), cscale_(1), nfval_(0), ngval_(0), ncval_(0),
      scaleLagrangian_(scaleLagrangian), HessianApprox_(HessianApprox) {

    fval_     = makePtr<ScalarController<Real,int>>();
    gradient_ = makePtr<VectorController<Real,int>>();
    conValue_ = makePtr<VectorController<Real,int>>();

    multiplier_    = dualConVec.clone();
    dualOptVector_ = dualOptVec.clone();
    dualConVector_ = dualConVec.clone();
    primConVector_ = primConVec.clone();
  }

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    obj_->update(x,type,iter);
    con_->update(x,type,iter);
    fval_->objectiveUpdate(type);
    gradient_->objectiveUpdate(type);
    conValue_->objectiveUpdate(type);
  }

  void setScaling(const Real fscale = 1.0, const Real cscale = 1.0) {
    fscale_ = fscale;
    cscale_ = cscale;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    // Compute objective function value
    Real val = getObjectiveValue(x,tol);
    val *= fscale_;
    // Compute penalty term
    const Real half(0.5);
    primConVector_->set(multiplier_->dual());
    primConVector_->axpy(half*cscale_*penaltyParameter_,*getConstraintVec(x,tol));
    val += cscale_*getConstraintVec(x,tol)->dot(*primConVector_);
    // Scale augmented Lagrangian
    if (scaleLagrangian_) {
      val /= penaltyParameter_;
    }
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Compute objective function gradient
    g.set(*getObjectiveGradient(x,tol));
    g.scale(fscale_);
    // Compute gradient of penalty
    dualConVector_->set(*multiplier_);
    dualConVector_->axpy(cscale_*penaltyParameter_,getConstraintVec(x,tol)->dual());
    con_->applyAdjointJacobian(*dualOptVector_,*dualConVector_,x,tol);
    g.axpy(cscale_,*dualOptVector_);
    // Compute gradient of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      const Real one(1);
      g.scale(one/penaltyParameter_);
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    hv.scale(fscale_);
    // Apply penalty Hessian to a vector
    if (HessianApprox_ < 3) {
      con_->applyJacobian(*primConVector_,v,x,tol);
      con_->applyAdjointJacobian(*dualOptVector_,primConVector_->dual(),x,tol);
      hv.axpy(cscale_*cscale_*penaltyParameter_,*dualOptVector_);
      if (HessianApprox_ == 1) {
        dualConVector_->set(*multiplier_);
        con_->applyAdjointHessian(*dualOptVector_,*dualConVector_,v,x,tol);
	hv.axpy(cscale_,*dualOptVector_);
      }
      if (HessianApprox_ == 0) {
        dualConVector_->set(*multiplier_);
        dualConVector_->axpy(cscale_*penaltyParameter_,getConstraintVec(x,tol)->dual());
        con_->applyAdjointHessian(*dualOptVector_,*dualConVector_,v,x,tol);
        hv.axpy(cscale_,*dualOptVector_);
      }
    }
    else {
      hv.zero();
    }
    // Build hessVec of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      hv.scale(static_cast<Real>(1)/penaltyParameter_);
    }
  }

  // Return objective function value
  Real getObjectiveValue(const Vector<Real> &x, Real &tol) {
    Real val(0);
    int key(0);
    bool isComputed = fval_->get(val,key);
    if ( !isComputed ) {
      val = obj_->value(x,tol); nfval_++;
      fval_->set(val,key);
    }
    return val;
  }

  // Compute objective function gradient
  const Ptr<const Vector<Real>> getObjectiveGradient(const Vector<Real> &x, Real &tol) {
    int key(0);
    if (!gradient_->isComputed(key)) {
      if (gradient_->isNull(key)) gradient_->allocate(*dualOptVector_,key);
      obj_->gradient(*gradient_->set(key),x,tol); ngval_++;
    }
    return gradient_->get(key);
  }

  // Return constraint value
  const Ptr<const Vector<Real>> getConstraintVec(const Vector<Real> &x, Real &tol) {
    int key(0);
    if (!conValue_->isComputed(key)) {
      if (conValue_->isNull(key)) conValue_->allocate(*primConVector_,key);
      con_->value(*conValue_->set(key),x,tol); ncval_++;
    }
    return conValue_->get(key);
  }

  // Return total number of constraint evaluations
  int getNumberConstraintEvaluations(void) const {
    return ncval_;
  }

  // Return total number of objective evaluations
  int getNumberFunctionEvaluations(void) const {
    return nfval_;
  }

  // Return total number of gradient evaluations
  int getNumberGradientEvaluations(void) const {
    return ngval_;
  }

  // Reset with upated penalty parameter
  void reset(const Vector<Real> &multiplier, const Real penaltyParameter) {
    nfval_ = 0; ngval_ = 0; ncval_ = 0;
    multiplier_->set(multiplier);
    penaltyParameter_ = penaltyParameter;
  }
}; // class AugmentedLagrangianObjective

} // namespace ROL

#endif
