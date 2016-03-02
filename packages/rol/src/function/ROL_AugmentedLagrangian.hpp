// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_AUGMENTEDLAGRANGIAN_H
#define ROL_AUGMENTEDLAGRANGIAN_H

#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::AugmentedLagrangian
    \brief Provides the interface to evaluate the augmented Lagrangian.

    ---
*/


namespace ROL {

template <class Real>
class AugmentedLagrangian : public Objective<Real> {
private:
  // Required for Augmented Lagrangian definition
  const Teuchos::RCP<Objective<Real> >          obj_;
  const Teuchos::RCP<EqualityConstraint<Real> > con_;
  Teuchos::RCP<Vector<Real> > lam_;
  Real penaltyParameter_;

  // Auxiliary storage
  Teuchos::RCP<Vector<Real> > dlam_;
  Teuchos::RCP<Vector<Real> > x_;
  Teuchos::RCP<Vector<Real> > dc1_;
  Teuchos::RCP<Vector<Real> > dc2_;

  // Objective and constraint evaluations
  Real fval_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > c_;

  // Evaluation counters
  int ncval_;
  int nfval_;
  int ngval_;

  // User defined options
  bool scaleLagrangian_;
  int HessianLevel_;

  // Flags to recompute quantities
  bool isConstraintComputed_;
  bool isValueComputed_;
  bool isGradientComputed_;

public:
  ~AugmentedLagrangian() {}

  AugmentedLagrangian(const Teuchos::RCP<Objective<Real> > &obj,
                      const Teuchos::RCP<EqualityConstraint<Real> > &con,
                      const Vector<Real> &lam,
                      const Real mu,
                      const Vector<Real> &x,
                      const Vector<Real> &c,
                      Teuchos::ParameterList &parlist)
    : obj_(obj), con_(con), penaltyParameter_(mu),
      fval_(0), ncval_(0), nfval_(0), ngval_(0),
      isConstraintComputed_(false), isValueComputed_(false), isGradientComputed_(false) {

    x_    = x.clone();
    g_    = x.dual().clone();
    dc1_  = x.dual().clone();
    dc2_  = c.clone();
    c_    = c.clone();
    lam_  = lam.clone();
    dlam_ = lam.clone();

    Teuchos::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_ = sublist.get("Use Scaled Augmented Lagrangian", false);
    HessianLevel_    = sublist.get("Level of Hessian Approximation",  0);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    con_->update(x,flag,iter);
    if ( flag ) {
      isConstraintComputed_ = false;
      isValueComputed_      = false;
      isGradientComputed_   = false;
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    if ( !isConstraintComputed_ ) {
      // Evaluate constraint
      con_->value(*c_,x,tol);
      ncval_++;
      isConstraintComputed_ = true;
    }
    if ( !isValueComputed_ ) {
      // Compute objective function value
      fval_ = obj_->value(x,tol);
      nfval_++;
      isValueComputed_ = true;
    }
    // Apply Lagrange multiplier to constraint
    Real cval = lam_->dot(c_->dual());
    // Compute penalty term
    Real pval = c_->dot(*c_);
    // Compute Augmented Lagrangian value
    Real val = 0.0, half = 0.5;
    if (scaleLagrangian_) {
      val = (fval_ + cval)/penaltyParameter_ + half*pval;
    }
    else {
      val = fval_ + cval + half*penaltyParameter_*pval;
    }
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    if ( !isConstraintComputed_ ) {
      // Evaluate constraint
      con_->value(*c_,x,tol);
      ncval_++;
      isConstraintComputed_ = true;
    }
    if ( !isGradientComputed_ ) {
      // Compute objective function gradient
      obj_->gradient(*g_,x,tol);
      ngval_++;
      isGradientComputed_ = true;
    }
    g.set(*g_);
    // Compute gradient of Augmented Lagrangian
    Real one = 1.0;
    dlam_->set(c_->dual());
    if ( scaleLagrangian_ ) {
      g.scale(one/penaltyParameter_);
      dlam_->axpy(one/penaltyParameter_,*lam_);
    }
    else {
      dlam_->scale(penaltyParameter_);
      dlam_->plus(*lam_);
    }
    con_->applyAdjointJacobian(*dc1_,*dlam_,x,tol);
    g.plus(*dc1_);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Real one = 1.0;
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    if (HessianLevel_ < 2) {
      con_->applyJacobian(*dc2_,v,x,tol);
      con_->applyAdjointJacobian(*dc1_,dc2_->dual(),x,tol);
      if (scaleLagrangian_) {
        hv.scale(one/penaltyParameter_);
        hv.plus(*dc1_);
      }
      else {
        hv.axpy(penaltyParameter_,*dc1_);
      }

      if (HessianLevel_ == 0) {
        if ( !isConstraintComputed_ ) {
          // Evaluate constraint
          con_->value(*c_,x,tol);
          ncval_++;
          isConstraintComputed_ = true;
        }
        // Apply Augmented Lagrangian Hessian to a vector
        dlam_->set(c_->dual());
        if ( scaleLagrangian_ ) {
          dlam_->axpy(one/penaltyParameter_,*lam_);
        }
        else {
          dlam_->scale(penaltyParameter_);
          dlam_->plus(*lam_);
        }
        con_->applyAdjointHessian(*dc1_,*dlam_,v,x,tol);
        hv.plus(*dc1_);
      }
    }
  }

  // Return objective function value
  Real getObjectiveValue(const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    if ( !isValueComputed_ ) {
      // Evaluate objective function value
      fval_ = obj_->value(x,tol);
      nfval_++;
      isValueComputed_ = true;
    }
    return fval_;
  }

  // Return objective function gradient
  void getObjectiveGradient(Vector<Real> &g, const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    if ( !isGradientComputed_ ) {
      // Compute objective function gradient
      obj_->gradient(*g_,x,tol);
      ngval_++;
      isGradientComputed_ = true;
    }
    g.set(*g_);
  }

  // Return constraint value
  void getConstraintVec(Vector<Real> &c, const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    if ( !isConstraintComputed_ ) {
      // Evaluate constraint
      con_->value(*c_,x,tol);
      ncval_++;
      isConstraintComputed_ = true;
    }
    c.set(*c_);
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
  void reset(const Vector<Real> &l, const Real penaltyParameter) {
    ncval_ = 0; nfval_ = 0; ngval_ = 0;
    lam_->set(l);
    penaltyParameter_ = penaltyParameter;
  }

}; // class AugmentedLagrangian

} // namespace ROL

#endif
