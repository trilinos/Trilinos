// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTEDLAGRANGIAN_SIMOPT_H
#define ROL_AUGMENTEDLAGRANGIAN_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_QuadraticPenalty_SimOpt.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::AugmentedLagrangian_SimOpt
    \brief Provides the interface to evaluate the SimOpt augmented Lagrangian.

    This class implements the SimOpt augmented Lagrangian functional for use with
    ROL::Reduced_AugmentedLagrangian_SimOpt.  Given a function
    \f$f:\mathcal{U}\times\mathcal{Z}\to\mathbb{R}\f$ and an equality constraint
    \f$c:\mathcal{U}\times\mathcal{Z}\to\mathcal{C}\f$, the augmented Lagrangian functional is
    \f[
       L_A(u,z,\lambda,\mu) = f(u,z) +
           \langle \lambda, c(u,z)\rangle_{\mathcal{C}^*,\mathcal{C}} +
           \frac{\mu}{2} \langle \mathfrak{R}c(u,z),c(u,z)\rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$\lambda\in\mathcal{C}^*\f$ denotes the Lagrange multiplier estimate,
    \f$\mu > 0\f$ is the penalty parameter and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ is the Riesz operator
    on the constraint space.

    This implementation permits the scaling of \f$L_A\f$ by \f$\mu^{-1}\f$ and also
    permits the Hessian approximation
    \f[
        \nabla^2_{uu} L_A(u,z,\lambda,\mu)v \approx \nabla^2_{uu} f(u,z) v
             + \mu c_u(u,z)^*\mathfrak{R} c_u(u,z)v,
        \quad
        \nabla^2_{uz} L_A(u,z,\lambda,\mu)v \approx \nabla^2_{uz} f(u,z) v
             + \mu c_u(u,z)^*\mathfrak{R} c_z(u,z)v,
    \f]
    \f[
        \nabla^2_{zu} L_A(u,z,\lambda,\mu)v \approx \nabla^2_{zu} f(u,z) v
             + \mu c_z(u,z)^*\mathfrak{R} c_u(u,z)v,
        \quad\text{and}\quad
        \nabla^2_{zz} L_A(u,z,\lambda,\mu)v \approx \nabla^2_{zz} f(u,z) v
             + \mu c_z(u,z)^*\mathfrak{R} c_z(u,z)v,
    \f]

    ---
*/


namespace ROL {

template <class Real>
class AugmentedLagrangian_SimOpt : public Objective_SimOpt<Real> {
private:
  // Required for Augmented Lagrangian definition
  const ROL::Ptr<Objective_SimOpt<Real> > obj_;
  ROL::Ptr<QuadraticPenalty_SimOpt<Real> > pen_;
  Real penaltyParameter_;

  // Auxiliary storage
  ROL::Ptr<Vector<Real> > dualSimVector_;
  ROL::Ptr<Vector<Real> > dualOptVector_;

  // Objective and constraint evaluations
  Real fval_;
  ROL::Ptr<Vector<Real> > gradient1_;
  ROL::Ptr<Vector<Real> > gradient2_;

  // Evaluation counters
  int nfval_;
  int ngval_;

  // User defined options
  bool scaleLagrangian_;

  // Flags to recompute quantities
  bool isValueComputed_;
  bool isGradient1Computed_;
  bool isGradient2Computed_;

public:
  AugmentedLagrangian_SimOpt(const ROL::Ptr<Objective_SimOpt<Real> > &obj,
                             const ROL::Ptr<Constraint_SimOpt<Real> > &con,
                             const Vector<Real> &multiplier,
                             const Real penaltyParameter,
                             const Vector<Real> &simVec,
                             const Vector<Real> &optVec,
                             const Vector<Real> &conVec,
                             ROL::ParameterList &parlist)
    : obj_(obj), penaltyParameter_(penaltyParameter),
      fval_(0), nfval_(0), ngval_(0), isValueComputed_(false),
      isGradient1Computed_(false), isGradient2Computed_(false) {

    gradient1_      = simVec.dual().clone();
    gradient2_      = optVec.dual().clone();
    dualSimVector_  = simVec.dual().clone();
    dualOptVector_  = optVec.dual().clone();

    ROL::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_  = sublist.get("Use Scaled Augmented Lagrangian", false);
    int HessianApprox = sublist.get("Level of Hessian Approximation",  0);

    pen_ = ROL::makePtr<QuadraticPenalty_SimOpt<Real>>(con,multiplier,penaltyParameter,simVec,optVec,conVec,scaleLagrangian_,HessianApprox);
  }

  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    obj_->update(u,z,flag,iter);
    pen_->update(u,z,flag,iter);
    isValueComputed_ = (flag ? false : isValueComputed_);
    isGradient1Computed_ = (flag ? false : isGradient1Computed_);
    isGradient2Computed_ = (flag ? false : isGradient2Computed_);
  }

  virtual Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    // Compute objective function value
    if ( !isValueComputed_ ) {
      fval_ = obj_->value(u,z,tol); nfval_++;
      isValueComputed_ = true;
    }
    // Compute penalty term
    Real pval = pen_->value(u,z,tol);
    // Compute augmented Lagrangian
    Real val = fval_;
    if (scaleLagrangian_) {
      val /= penaltyParameter_;
    }
    return val + pval;
  }

  virtual void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    // Compute objective function gradient
    if ( !isGradient1Computed_ ) {
      obj_->gradient_1(*gradient1_,u,z,tol); ngval_++;
      isGradient1Computed_ = true;
    }
    g.set(*gradient1_);
    // Compute gradient of penalty
    pen_->gradient_1(*dualSimVector_,u,z,tol);
    // Compute gradient of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      g.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    g.plus(*dualSimVector_);
  }

  virtual void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    // Compute objective function gradient
    if ( !isGradient2Computed_ ) {
      obj_->gradient_2(*gradient2_,u,z,tol); ngval_++;
      isGradient2Computed_ = true;
    }
    g.set(*gradient2_);
    // Compute gradient of penalty
    pen_->gradient_2(*dualOptVector_,u,z,tol);
    // Compute gradient of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      g.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    g.plus(*dualOptVector_);
  }

  virtual void hessVec_11( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec_11(hv,v,u,z,tol);
    // Apply penalty Hessian to a vector
    pen_->hessVec_11(*dualSimVector_,v,u,z,tol);
    // Build hessVec of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      hv.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    hv.plus(*dualSimVector_);
  }

  virtual void hessVec_12( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec_12(hv,v,u,z,tol);
    // Apply penalty Hessian to a vector
    pen_->hessVec_12(*dualSimVector_,v,u,z,tol);
    // Build hessVec of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      hv.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    hv.plus(*dualSimVector_);
  }

  virtual void hessVec_21( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec_21(hv,v,u,z,tol);
    // Apply penalty Hessian to a vector
    pen_->hessVec_21(*dualOptVector_,v,u,z,tol);
    // Build hessVec of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      hv.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    hv.plus(*dualOptVector_);
  }

  virtual void hessVec_22( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec_22(hv,v,u,z,tol);
    // Apply penalty Hessian to a vector
    pen_->hessVec_22(*dualOptVector_,v,u,z,tol);
    // Build hessVec of Augmented Lagrangian
    if ( scaleLagrangian_ ) {
      hv.scale(static_cast<Real>(1)/penaltyParameter_);
    }
    hv.plus(*dualOptVector_);
  }

  // Return objective function value
  virtual Real getObjectiveValue(const Vector<Real> &u, const Vector<Real> &z) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Evaluate objective function value
    if ( !isValueComputed_ ) {
      fval_ = obj_->value(u,z,tol); nfval_++;
      isValueComputed_ = true;
    }
    return fval_;
  }

  // Return constraint value
  virtual void getConstraintVec(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z) {
    pen_->getConstraintVec(c,u,z);
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
