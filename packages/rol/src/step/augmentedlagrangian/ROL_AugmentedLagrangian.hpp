// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTEDLAGRANGIAN_H
#define ROL_AUGMENTEDLAGRANGIAN_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_QuadraticPenalty.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::AugmentedLagrangian
    \brief Provides the interface to evaluate the augmented Lagrangian.

    This class implements the augmented Lagrangian functional for use with
    ROL::AugmentedLagrangianStep.  Given a function
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
class AugmentedLagrangian : public Objective<Real> {
private:
  // Required for Augmented Lagrangian definition
  const ROL::Ptr<Objective<Real> > obj_;
  ROL::Ptr<QuadraticPenalty<Real> > pen_;
  Real penaltyParameter_;

  // Auxiliary storage
  ROL::Ptr<Vector<Real> > dualOptVector_;

  // Objective and constraint evaluations
  Real fval_;
  ROL::Ptr<Vector<Real> > gradient_;

  // Objective function scaling
  Real fscale_;

  // Evaluation counters
  int nfval_;
  int ngval_;

  // User defined options
  bool scaleLagrangian_;

  // Flags to recompute quantities
  bool isValueComputed_;
  bool isGradientComputed_;

public:
  /** \brief Constructor.

      This creates a valid AugmentedLagrangian object.
      @param[in]          obj              is an objective function.
      @param[in]          con              is an equality constraint.
      @param[in]          mulitplier       is a Lagrange multiplier vector.
      @param[in]          penaltyParameter is the penalty parameter.
      @param[in]          optVec           is an optimization space vector.
      @param[in]          conVec           is a constraint space vector.
      @param[in]          parlist          is a parameter list.
  */
  AugmentedLagrangian(const ROL::Ptr<Objective<Real> > &obj,
                      const ROL::Ptr<Constraint<Real> > &con,
                      const Vector<Real> &multiplier,
                      const Real penaltyParameter,
                      const Vector<Real> &optVec,
                      const Vector<Real> &conVec,
                      ROL::ParameterList &parlist)
    : obj_(obj), penaltyParameter_(penaltyParameter),
      fval_(0), fscale_(1),
      nfval_(0), ngval_(0),
      isValueComputed_(false), isGradientComputed_(false) {

    gradient_      = optVec.dual().clone();
    dualOptVector_ = optVec.dual().clone();

    ROL::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_  = sublist.get("Use Scaled Augmented Lagrangian", false);
    int HessianApprox = sublist.get("Level of Hessian Approximation",  0);

    pen_ = ROL::makePtr<QuadraticPenalty<Real>>(con,multiplier,penaltyParameter,optVec,conVec,scaleLagrangian_,HessianApprox);
  }

  /** \brief Constructor.

      This creates a valid AugmentedLagrangian object.
      @param[in]          obj              is an objective function.
      @param[in]          con              is an equality constraint.
      @param[in]          mulitplier       is a Lagrange multiplier vector.
      @param[in]          penaltyParameter is the penalty parameter.
      @param[in]          optVec           is an optimization space vector.
      @param[in]          conVec           is a constraint space vector.
      @param[in]          parlist          is a parameter list.
  */
  AugmentedLagrangian(const ROL::Ptr<Objective<Real> > &obj,
                      const ROL::Ptr<Constraint<Real> > &con,
                      const Vector<Real> &multiplier,
                      const Real penaltyParameter,
                      const Vector<Real> &optVec,
                      const Vector<Real> &conVec,
                      const bool scaleLagrangian,
                      const int HessianApprox)
    : obj_(obj), penaltyParameter_(penaltyParameter),
      fval_(0), fscale_(1),
      nfval_(0), ngval_(0),
      scaleLagrangian_(scaleLagrangian),
      isValueComputed_(false), isGradientComputed_(false) {

    gradient_      = optVec.dual().clone();
    dualOptVector_ = optVec.dual().clone();

    pen_ = ROL::makePtr<QuadraticPenalty<Real>>(con,multiplier,penaltyParameter,optVec,conVec,scaleLagrangian_,HessianApprox);
  }

  /** \brief Null constructor.

      This constructor is only used for inheritance and does not create a
      valid AugmentedLagrangian object.  Do not use.
  */
  AugmentedLagrangian()
   : obj_(ROL::nullPtr), pen_(ROL::nullPtr), dualOptVector_(ROL::nullPtr),
     fval_(0), gradient_(ROL::nullPtr), fscale_(1),
     nfval_(0), ngval_(0),
     scaleLagrangian_(false), isValueComputed_(false), isGradientComputed_(false) {}

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    pen_->update(x,flag,iter);
    isValueComputed_ = ((flag || (!flag && iter < 0)) ? false : isValueComputed_);
    isGradientComputed_ = ((flag || (!flag && iter < 0)) ? false : isGradientComputed_);
  }

  void setScaling(const Real fscale, const Real cscale = 1.0) {
    fscale_ = fscale;
    pen_->setScaling(cscale);
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
    // Compute objective function value
    if ( !isValueComputed_ ) {
      fval_ = obj_->value(x,tol); nfval_++;
      isValueComputed_ = true;
    }
    // Compute penalty term
    Real pval = pen_->value(x,tol);
    // Compute augmented Lagrangian
    Real val = fscale_*fval_;
    if (scaleLagrangian_) {
      val /= penaltyParameter_;
    }
    return val + pval;
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Compute objective function gradient
    if ( !isGradientComputed_ ) {
      obj_->gradient(*gradient_,x,tol); ngval_++;
      isGradientComputed_ = true;
    }
    g.set(*gradient_);
    g.scale(fscale_);
    // Compute gradient of penalty
    pen_->gradient(*dualOptVector_,x,tol);
    // Compute gradient of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      g.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    g.plus(*dualOptVector_);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    hv.scale(fscale_);
    // Apply penalty Hessian to a vector
    pen_->hessVec(*dualOptVector_,v,x,tol);
    // Build hessVec of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      hv.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    hv.plus(*dualOptVector_);
  }

  // Return objective function value
  virtual Real getObjectiveValue(const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Evaluate objective function value
    if ( !isValueComputed_ ) {
      fval_ = obj_->value(x,tol); nfval_++;
      isValueComputed_ = true;
    }
    return fval_;
  }

  const Ptr<const Vector<Real>> getObjectiveGradient(const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Compute objective function gradient
    if ( !isGradientComputed_ ) {
      obj_->gradient(*gradient_,x,tol); ngval_++;
      isGradientComputed_ = true;
    }
    return gradient_;
  }

  // Return constraint value
  virtual void getConstraintVec(Vector<Real> &c, const Vector<Real> &x) {
    pen_->getConstraintVec(c,x);
  }

  // Return total number of constraint evaluations
  virtual int getNumberConstraintEvaluations(void) const {
    return pen_->getNumberConstraintEvaluations();
  }

  // Return total number of objective evaluations
  virtual int getNumberFunctionEvaluations(void) const {
    return nfval_;
  }

  // Return total number of gradient evaluations
  virtual int getNumberGradientEvaluations(void) const {
    return ngval_;
  }

  // Reset with upated penalty parameter
  virtual void reset(const Vector<Real> &multiplier, const Real penaltyParameter) {
    nfval_ = 0; ngval_ = 0;
    pen_->reset(multiplier,penaltyParameter);
  }
}; // class AugmentedLagrangian

} // namespace ROL

#endif
