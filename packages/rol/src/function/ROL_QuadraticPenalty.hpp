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

#ifndef ROL_QUADRATICPENALTY_H
#define ROL_QUADRATICPENALTY_H

#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::QuadraticPenalty
    \brief Provides the interface to evaluate the quadratic constraint penalty.

    ---
*/


namespace ROL {

template <class Real>
class QuadraticPenalty : public Objective<Real> {
private:
  // Required for quadratic penalty definition
  const Teuchos::RCP<EqualityConstraint<Real> > con_;
  Teuchos::RCP<Vector<Real> > multiplier_;
  Real penaltyParameter_;

  // Auxiliary storage
  Teuchos::RCP<Vector<Real> > primalMultiplierVector_;
  Teuchos::RCP<Vector<Real> > dualOptVector_;
  Teuchos::RCP<Vector<Real> > primalConVector_;

  // Constraint evaluations
  Teuchos::RCP<Vector<Real> > conValue_;

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
  QuadraticPenalty(const Teuchos::RCP<EqualityConstraint<Real> > &con,
                   const Vector<Real> &multiplier,
                   const Real penaltyParameter,
                   const Vector<Real> &optVec,
                   const Vector<Real> &conVec,
                   const bool useScaling = false,
                   const int HessianApprox = 0 )
    : con_(con), penaltyParameter_(penaltyParameter), ncval_(0),
      useScaling_(useScaling), HessianApprox_(HessianApprox), isConstraintComputed_(false) {

    dualOptVector_          = optVec.dual().clone();
    primalConVector_        = conVec.clone();
    conValue_               = conVec.clone();
    multiplier_             = multiplier.clone();
    primalMultiplierVector_ = multiplier.clone();
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
    isConstraintComputed_ = ( flag ? false : isConstraintComputed_ );
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(x,tol);
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

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Evaluate constraint
    evaluateConstraint(x,tol);
    // Compute gradient of Augmented Lagrangian
    const Real one(1);
    primalMultiplierVector_->set(conValue_->dual());
    if ( useScaling_ ) {
      primalMultiplierVector_->axpy(one/penaltyParameter_,*multiplier_);
    }
    else {
      primalMultiplierVector_->scale(penaltyParameter_);
      primalMultiplierVector_->plus(*multiplier_);
    }
    con_->applyAdjointJacobian(g,*primalMultiplierVector_,x,tol);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    if (HessianApprox_ < 2) {
      con_->applyJacobian(*primalConVector_,v,x,tol);
      con_->applyAdjointJacobian(hv,primalConVector_->dual(),x,tol);
      if (!useScaling_) {
        hv.scale(penaltyParameter_);
      }

      if (HessianApprox_ == 0) {
        // Evaluate constraint
        evaluateConstraint(x,tol);
        // Apply Augmented Lagrangian Hessian to a vector
        const Real one(1);
        primalMultiplierVector_->set(conValue_->dual());
        if ( useScaling_ ) {
          primalMultiplierVector_->axpy(one/penaltyParameter_,*multiplier_);
        }
        else {
          primalMultiplierVector_->scale(penaltyParameter_);
          primalMultiplierVector_->plus(*multiplier_);
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
