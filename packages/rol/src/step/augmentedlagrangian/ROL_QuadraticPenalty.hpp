// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATICPENALTY_H
#define ROL_QUADRATICPENALTY_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::QuadraticPenalty
    \brief Provides the interface to evaluate the quadratic constraint penalty.

    This class implements the quadratic constraint penalty functional.
    Given an equality constraint \f$c:\mathcal{X}\to\mathcal{C}\f$, the
    quadratic penalty functional is
    \f[
       Q(x,\lambda,\mu) =
           \langle \lambda, c(x)\rangle_{\mathcal{C}^*,\mathcal{C}} +
           \frac{\mu}{2} \langle \mathfrak{R}c(x),c(x)\rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$\lambda\in\mathcal{C}^*\f$ denotes a multiplier,
    \f$\mu > 0\f$ is the penalty parameter and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ is the Riesz operator
    on the constraint space.

    This implementation permits the scaling of \f$Q\f$ by \f$\mu^{-1}\f$ and also
    permits the Hessian approximation
    \f[
        \nabla^2_x Q(x,\lambda,\mu)v \approx \mu c'(x)^*\mathfrak{R} c'(x)v.
    \f]

    ---
*/


namespace ROL {

template <class Real>
class QuadraticPenalty : public Objective<Real> {
private:
  // Required for quadratic penalty definition
  const ROL::Ptr<Constraint<Real> > con_;
  ROL::Ptr<Vector<Real> > multiplier_;
  Real penaltyParameter_;

  // Auxiliary storage
  ROL::Ptr<Vector<Real> > primalMultiplierVector_;
  ROL::Ptr<Vector<Real> > dualOptVector_;
  ROL::Ptr<Vector<Real> > primalConVector_;

  // Constraint evaluations
  ROL::Ptr<Vector<Real> > conValue_;
  Real cscale_;

  // Evaluation counters
  int ncval_;

  // User defined options
  const bool useScaling_;
  const int HessianApprox_;

  // Flags to recompute quantities
  bool isConstraintComputed_;

  void evaluateConstraint(const Vector<Real> &x, Real &tol) {
    if ( !isConstraintComputed_ ) {
      // Evaluate constraint
      con_->value(*conValue_,x,tol); ncval_++;
      isConstraintComputed_ = true;
    }
  }

public:
  QuadraticPenalty(const ROL::Ptr<Constraint<Real> > &con,
                   const Vector<Real> &multiplier,
                   const Real penaltyParameter,
                   const Vector<Real> &optVec,
                   const Vector<Real> &conVec,
                   const bool useScaling = false,
                   const int HessianApprox = 0 )
    : con_(con), penaltyParameter_(penaltyParameter), cscale_(1), ncval_(0),
      useScaling_(useScaling), HessianApprox_(HessianApprox), isConstraintComputed_(false) {

    dualOptVector_          = optVec.dual().clone();
    primalConVector_        = conVec.clone();
    conValue_               = conVec.clone();
    multiplier_             = multiplier.clone();
    primalMultiplierVector_ = multiplier.clone();
  }

  void setScaling(const Real cscale = 1) {
    cscale_ = cscale;
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
    isConstraintComputed_ = ((flag || (!flag && iter < 0)) ? false : isConstraintComputed_);
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(x,tol);
    // Apply Lagrange multiplier to constraint
    Real cval = cscale_*multiplier_->dot(conValue_->dual());
    // Compute penalty term
    Real pval = cscale_*cscale_*conValue_->dot(*conValue_);
    // Compute quadratic penalty value
    const Real half(0.5);
    Real val(0);
    if (useScaling_) {
      val = cval/penaltyParameter_ + half*pval;
    }
    else {
      val = cval + half*penaltyParameter_*pval;
    }
    return val;
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(x,tol);
    // Compute gradient of Augmented Lagrangian
    primalMultiplierVector_->set(conValue_->dual());
    if ( useScaling_ ) {
      primalMultiplierVector_->scale(cscale_*cscale_);
      primalMultiplierVector_->axpy(cscale_/penaltyParameter_,*multiplier_);
    }
    else {
      primalMultiplierVector_->scale(cscale_*cscale_*penaltyParameter_);
      primalMultiplierVector_->axpy(cscale_,*multiplier_);
    }
    con_->applyAdjointJacobian(g,*primalMultiplierVector_,x,tol);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    if (HessianApprox_ < 3) {
      con_->update(x);
      con_->applyJacobian(*primalConVector_,v,x,tol);
      con_->applyAdjointJacobian(hv,primalConVector_->dual(),x,tol);
      if (!useScaling_) {
        hv.scale(cscale_*cscale_*penaltyParameter_);
      }
      else {
        hv.scale(cscale_*cscale_);
      }

      if (HessianApprox_ == 1) {
        // Apply Augmented Lagrangian Hessian to a vector
        primalMultiplierVector_->set(*multiplier_);
        if ( useScaling_ ) {
          primalMultiplierVector_->scale(cscale_/penaltyParameter_);
        }
        else {
          primalMultiplierVector_->scale(cscale_);
        }
        con_->applyAdjointHessian(*dualOptVector_,*primalMultiplierVector_,v,x,tol);
        hv.plus(*dualOptVector_);
      }

      if (HessianApprox_ == 0) {
        // Evaluate constraint
        evaluateConstraint(x,tol);
        // Apply Augmented Lagrangian Hessian to a vector
        primalMultiplierVector_->set(conValue_->dual());
        if ( useScaling_ ) {
          primalMultiplierVector_->scale(cscale_*cscale_);
          primalMultiplierVector_->axpy(cscale_/penaltyParameter_,*multiplier_);
        }
        else {
          primalMultiplierVector_->scale(cscale_*cscale_*penaltyParameter_);
          primalMultiplierVector_->axpy(cscale_,*multiplier_);
        }
        con_->applyAdjointHessian(*dualOptVector_,*primalMultiplierVector_,v,x,tol);
        hv.plus(*dualOptVector_);
      }
    }
    else {
      hv.zero();
    }
  }

  // Return constraint value
  virtual void getConstraintVec(Vector<Real> &c, const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Evaluate constraint
    evaluateConstraint(x,tol);
    c.set(*conValue_);
  }

  // Return total number of constraint evaluations
  virtual int getNumberConstraintEvaluations(void) const {
    return ncval_;
  }

  // Reset with upated penalty parameter
  virtual void reset(const Vector<Real> &multiplier, const Real penaltyParameter) {
    ncval_ = 0;
    multiplier_->set(multiplier);
    penaltyParameter_ = penaltyParameter;
  }
}; // class AugmentedLagrangian

} // namespace ROL

#endif
