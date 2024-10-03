// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATICPENALTY_SIMOPT_H
#define ROL_QUADRATICPENALTY_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::QuadraticPenalty_SimOpt
    \brief Provides the interface to evaluate the quadratic SimOpt constraint penalty.

    This class implements the quadratic SimOpt constraint penalty functional.
    Given an equality constraint \f$c:\mathcal{U}\times\mathcal{Z}\to\mathcal{C}\f$, the
    quadratic penalty functional is
    \f[
       Q(u,z,\lambda,\mu) =
           \langle \lambda, c(u,z)\rangle_{\mathcal{C}^*,\mathcal{C}} +
           \frac{\mu}{2} \langle \mathfrak{R}c(u,z),c(u,z)\rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$\lambda\in\mathcal{C}^*\f$ denotes a multiplier,
    \f$\mu > 0\f$ is the penalty parameter and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ is the Riesz operator
    on the constraint space.

    This implementation permits the scaling of \f$Q\f$ by \f$\mu^{-1}\f$ and also
    permits the Hessian approximation
    \f[
        \nabla^2_{uu} Q(u,z,\lambda,\mu)v \approx \mu c_u(u,z)^*\mathfrak{R} c_u(u,z)v,
        \quad
        \nabla^2_{uz} Q(u,z,\lambda,\mu)v \approx \mu c_u(u,z)^*\mathfrak{R} c_z(u,z)v,
    \f]
    \f[
        \nabla^2_{zu} Q(u,z,\lambda,\mu)v \approx \mu c_z(u,z)^*\mathfrak{R} c_u(u,z)v,
        \quad\text{and}\quad
        \nabla^2_{zz} Q(u,z,\lambda,\mu)v \approx \mu c_z(u,z)^*\mathfrak{R} c_z(u,z)v.
    \f]

    ---
*/


namespace ROL {

template <class Real>
class QuadraticPenalty_SimOpt : public Objective_SimOpt<Real> {
private:
  // Required for quadratic penalty definition
  const ROL::Ptr<Constraint_SimOpt<Real> > con_;
  ROL::Ptr<Vector<Real> > multiplier_;
  Real penaltyParameter_;

  // Auxiliary storage
  ROL::Ptr<Vector<Real> > primalMultiplierVector_;
  ROL::Ptr<Vector<Real> > dualSimVector_;
  ROL::Ptr<Vector<Real> > dualOptVector_;
  ROL::Ptr<Vector<Real> > primalConVector_;

  // Constraint evaluations
  ROL::Ptr<Vector<Real> > conValue_;

  // Evaluation counters
  int ncval_;

  // User defined options
  const bool useScaling_;
  const int HessianApprox_;

  // Flags to recompute quantities
  bool isConstraintComputed_;

  void evaluateConstraint(const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    if ( !isConstraintComputed_ ) {
      // Evaluate constraint
      con_->value(*conValue_,u,z,tol); ncval_++;
      isConstraintComputed_ = true;
    }
  }

public:
  QuadraticPenalty_SimOpt(const ROL::Ptr<Constraint_SimOpt<Real> > &con,
                          const Vector<Real> &multiplier,
                          const Real penaltyParameter,
                          const Vector<Real> &simVec,
                          const Vector<Real> &optVec,
                          const Vector<Real> &conVec,
                          const bool useScaling = false,
                          const int HessianApprox = 0 )
    : con_(con), penaltyParameter_(penaltyParameter), ncval_(0),
      useScaling_(useScaling), HessianApprox_(HessianApprox), isConstraintComputed_(false) {

    dualSimVector_          = simVec.dual().clone();
    dualOptVector_          = optVec.dual().clone();
    primalConVector_        = conVec.clone();
    conValue_               = conVec.clone();
    multiplier_             = multiplier.clone();
    primalMultiplierVector_ = multiplier.clone();
  }

  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    con_->update(u,z,flag,iter);
    isConstraintComputed_ = ( flag ? false : isConstraintComputed_ );
  }

  virtual Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(u,z,tol);
    // Apply Lagrange multiplier to constraint
    Real cval = multiplier_->dot(conValue_->dual());
    // Compute penalty term
    Real pval = conValue_->dot(*conValue_);
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

  virtual void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(u,z,tol);
    // Compute gradient of Augmented Lagrangian
    primalMultiplierVector_->set(conValue_->dual());
    if ( useScaling_ ) {
      primalMultiplierVector_->axpy(static_cast<Real>(1)/penaltyParameter_,*multiplier_);
    }
    else {
      primalMultiplierVector_->scale(penaltyParameter_);
      primalMultiplierVector_->plus(*multiplier_);
    }
    con_->applyAdjointJacobian_1(g,*primalMultiplierVector_,u,z,tol);
  }

  virtual void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(u,z,tol);
    // Compute gradient of Augmented Lagrangian
    primalMultiplierVector_->set(conValue_->dual());
    if ( useScaling_ ) {
      primalMultiplierVector_->axpy(static_cast<Real>(1)/penaltyParameter_,*multiplier_);
    }
    else {
      primalMultiplierVector_->scale(penaltyParameter_);
      primalMultiplierVector_->plus(*multiplier_);
    }
    con_->applyAdjointJacobian_2(g,*primalMultiplierVector_,u,z,tol);
  }

  virtual void hessVec_11( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    if (HessianApprox_ < 2) {
      con_->applyJacobian_1(*primalConVector_,v,u,z,tol);
      con_->applyAdjointJacobian_1(hv,primalConVector_->dual(),u,z,tol);
      if (!useScaling_) {
        hv.scale(penaltyParameter_);
      }

      if (HessianApprox_ == 0) {
        // Evaluate constraint
        evaluateConstraint(u,z,tol);
        // Apply Augmented Lagrangian Hessian to a vector
        primalMultiplierVector_->set(conValue_->dual());
        if ( useScaling_ ) {
          primalMultiplierVector_->axpy(static_cast<Real>(1)/penaltyParameter_,*multiplier_);
        }
        else {
          primalMultiplierVector_->scale(penaltyParameter_);
          primalMultiplierVector_->plus(*multiplier_);
        }
        con_->applyAdjointHessian_11(*dualSimVector_,*primalMultiplierVector_,v,u,z,tol);
        hv.plus(*dualSimVector_);
      }
    }
    else {
      hv.zero();
    }
  }

  virtual void hessVec_12( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    if (HessianApprox_ < 2) {
      con_->applyJacobian_2(*primalConVector_,v,u,z,tol);
      con_->applyAdjointJacobian_1(hv,primalConVector_->dual(),u,z,tol);
      if (!useScaling_) {
        hv.scale(penaltyParameter_);
      }

      if (HessianApprox_ == 0) {
        // Evaluate constraint
        evaluateConstraint(u,z,tol);
        // Apply Augmented Lagrangian Hessian to a vector
        primalMultiplierVector_->set(conValue_->dual());
        if ( useScaling_ ) {
          primalMultiplierVector_->axpy(static_cast<Real>(1)/penaltyParameter_,*multiplier_);
        }
        else {
          primalMultiplierVector_->scale(penaltyParameter_);
          primalMultiplierVector_->plus(*multiplier_);
        }
        con_->applyAdjointHessian_21(*dualSimVector_,*primalMultiplierVector_,v,u,z,tol);
        hv.plus(*dualSimVector_);
      }
    }
    else {
      hv.zero();
    }
  }

  virtual void hessVec_21( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    if (HessianApprox_ < 2) {
      con_->applyJacobian_1(*primalConVector_,v,u,z,tol);
      con_->applyAdjointJacobian_2(hv,primalConVector_->dual(),u,z,tol);
      if (!useScaling_) {
        hv.scale(penaltyParameter_);
      }

      if (HessianApprox_ == 0) {
        // Evaluate constraint
        evaluateConstraint(u,z,tol);
        // Apply Augmented Lagrangian Hessian to a vector
        primalMultiplierVector_->set(conValue_->dual());
        if ( useScaling_ ) {
          primalMultiplierVector_->axpy(static_cast<Real>(1)/penaltyParameter_,*multiplier_);
        }
        else {
          primalMultiplierVector_->scale(penaltyParameter_);
          primalMultiplierVector_->plus(*multiplier_);
        }
        con_->applyAdjointHessian_12(*dualOptVector_,*primalMultiplierVector_,v,u,z,tol);
        hv.plus(*dualOptVector_);
      }
    }
    else {
      hv.zero();
    }
  }

  virtual void hessVec_22( Vector<Real> &hv, const Vector<Real> &v,
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    // Apply objective Hessian to a vector
    if (HessianApprox_ < 2) {
      con_->applyJacobian_2(*primalConVector_,v,u,z,tol);
      con_->applyAdjointJacobian_2(hv,primalConVector_->dual(),u,z,tol);
      if (!useScaling_) {
        hv.scale(penaltyParameter_);
      }

      if (HessianApprox_ == 0) {
        // Evaluate constraint
        evaluateConstraint(u,z,tol);
        // Apply Augmented Lagrangian Hessian to a vector
        primalMultiplierVector_->set(conValue_->dual());
        if ( useScaling_ ) {
          primalMultiplierVector_->axpy(static_cast<Real>(1)/penaltyParameter_,*multiplier_);
        }
        else {
          primalMultiplierVector_->scale(penaltyParameter_);
          primalMultiplierVector_->plus(*multiplier_);
        }
        con_->applyAdjointHessian_22(*dualOptVector_,*primalMultiplierVector_,v,u,z,tol);
        hv.plus(*dualOptVector_);
      }
    }
    else {
      hv.zero();
    }
  }

  // Return constraint value
  virtual void getConstraintVec(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Evaluate constraint
    evaluateConstraint(u,z,tol);
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
