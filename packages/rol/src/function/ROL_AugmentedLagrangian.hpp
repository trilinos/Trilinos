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
#include "ROL_QuadraticPenalty.hpp"
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
  const Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<QuadraticPenalty<Real> > pen_;
  Real penaltyParameter_;

  // Auxiliary storage
  Teuchos::RCP<Vector<Real> > dualOptVector_;

  // Objective and constraint evaluations
  Real fval_;
  Teuchos::RCP<Vector<Real> > gradient_;

  // Evaluation counters
  int nfval_;
  int ngval_;

  // User defined options
  bool scaleLagrangian_;

  // Flags to recompute quantities
  bool isValueComputed_;
  bool isGradientComputed_;

public:
  AugmentedLagrangian(const Teuchos::RCP<Objective<Real> > &obj,
                      const Teuchos::RCP<EqualityConstraint<Real> > &con,
                      const Vector<Real> &multiplier,
                      const Real penaltyParameter,
                      const Vector<Real> &optVec,
                      const Vector<Real> &conVec,
                      Teuchos::ParameterList &parlist)
    : obj_(obj), penaltyParameter_(penaltyParameter),
      fval_(0), nfval_(0), ngval_(0), isValueComputed_(false), isGradientComputed_(false) {

    gradient_      = optVec.dual().clone();
    dualOptVector_ = optVec.dual().clone();

    Teuchos::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_  = sublist.get("Use Scaled Augmented Lagrangian", false);
    int HessianApprox = sublist.get("Level of Hessian Approximation",  0);

    pen_ = Teuchos::rcp(new QuadraticPenalty<Real>(con,multiplier,penaltyParameter,optVec,conVec,scaleLagrangian_,HessianApprox));
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    pen_->update(x,flag,iter);
    isValueComputed_ = (flag ? false : isValueComputed_);
    isGradientComputed_ = (flag ? false : isGradientComputed_);
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
    Real val = fval_;
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

  // Return objective function gradient
  virtual void getObjectiveGradient(Vector<Real> &g, const Vector<Real> &x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Compute objective function gradient
    if ( !isGradientComputed_ ) {
      obj_->gradient(*gradient_,x,tol); ngval_++;
      isGradientComputed_ = true;
    }
    g.set(*gradient_);
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
