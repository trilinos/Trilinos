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
  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<EqualityConstraint<Real> > con_;
  Teuchos::RCP<Vector<Real> > lam_;
  Teuchos::RCP<Vector<Real> > dlam_;
  Teuchos::RCP<Vector<Real> > x_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > c_;
  Teuchos::RCP<Vector<Real> > dc1_;
  Teuchos::RCP<Vector<Real> > dc2_;

  Real fval_;
  int ncval_;
  int nfval_;
  int ngval_;

  Real penaltyParameter_;
  bool scaleLagrangian_;
  int HessianLevel_;

  bool isConstraintComputed_;
  bool isValueComputed_;
  bool isGradientComputed_;

public:
  ~AugmentedLagrangian() {}

  AugmentedLagrangian(Objective<Real> &obj, EqualityConstraint<Real> &con, 
                const Vector<Real> &x, const Vector<Real> &c,
                const Vector<Real> &l, const Real mu,
                Teuchos::ParameterList &parlist)
    : fval_(0.0), ncval_(0), nfval_(0), ngval_(0), penaltyParameter_(mu),
      isConstraintComputed_(false), isValueComputed_(false), isGradientComputed_(false) {
    obj_ = Teuchos::rcp(&obj, false);
    con_ = Teuchos::rcp(&con, false);

    x_    = x.clone();
    g_    = x.dual().clone();
    dc1_  = x.dual().clone();
    dc2_  = c.clone();
    c_    = c.clone();
    lam_  = l.clone();
    dlam_ = l.clone();
    lam_->set(l);

    Teuchos::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_ = sublist.get("Use Scaled Augmented Lagrangian", false);
    HessianLevel_    = sublist.get("Level of Hessian Approximation",  0);
  }

  Real getObjectiveValue(const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON);
    if ( !isValueComputed_ ) {
      // Evaluate objective function value
      fval_ = obj_->value(x,tol);
      nfval_++;
      isValueComputed_ = true;
    }
    return fval_;
  }

  void getObjectiveGradient(Vector<Real> &g, const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON);
    if ( !isGradientComputed_ ) {
      // Compute objective function gradient
      obj_->gradient(*g_,x,tol);
      ngval_++;
      isGradientComputed_ = true;
    }
    g.set(*g_);
  }

  void getConstraintVec(Vector<Real> &c, const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON);
    if ( !isConstraintComputed_ ) {
      // Evaluate constraint
      con_->value(*c_,x,tol);
      ncval_++;
      isConstraintComputed_ = true;
    }
    c.set(*c_);
  }

  int getNumberConstraintEvaluations(void) {
    return ncval_;
  }

  int getNumberFunctionEvaluations(void) {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) {
    return ngval_;
  }

  void reset(const Vector<Real> &lam, const Real penaltyParameter) {
    ncval_ = 0.; nfval_ = 0.; ngval_ = 0.;
    lam_->set(lam);
    penaltyParameter_ = penaltyParameter;
  }

  /** \brief Update augmented Lagrangian function. 

      This function updates the augmented Lagrangian function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    con_->update(x,flag,iter);
    if ( flag ) {
      isConstraintComputed_ = false;
      isValueComputed_      = false;
      isGradientComputed_   = false;
    }
  }

  /** \brief Compute value.

      This function returns the augmented Lagrangian value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact augmented Lagrangian computation.
  */
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
    Real val = 0.0;
    if (scaleLagrangian_) {
      val = (fval_ + cval)/penaltyParameter_ + 0.5*pval;
    }
    else {
      val = fval_ + cval + 0.5*penaltyParameter_*pval;
    }
    return val;
  }

  /** \brief Compute gradient.

      This function returns the augmented Lagrangian gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact augmented Lagrangian computation.
  */
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
    dlam_->set(c_->dual());
    if ( scaleLagrangian_ ) {
      g.scale(1./penaltyParameter_);
      dlam_->axpy(1./penaltyParameter_,*lam_);
    }
    else {
      dlam_->scale(penaltyParameter_);
      dlam_->plus(*lam_);
    }
    con_->applyAdjointJacobian(*dc1_,*dlam_,x,tol);
    g.plus(*dc1_);
  }

  /** \brief Apply Hessian approximation to vector.

      This function applies the Hessian of the augmented Lagrangian to the vector \f$v\f$.
      @param[out]         hv  is the the action of the Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact augmented Lagrangian computation.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    if (HessianLevel_ < 2) {
      con_->applyJacobian(*dc2_,v,x,tol);
      con_->applyAdjointJacobian(*dc1_,dc2_->dual(),x,tol);
      if (scaleLagrangian_) {
        hv.scale(1./penaltyParameter_);
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
          dlam_->axpy(1./penaltyParameter_,*lam_);
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

}; // class AugmentedLagrangian

} // namespace ROL

#endif
