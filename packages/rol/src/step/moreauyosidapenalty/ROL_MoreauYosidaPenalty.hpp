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

#ifndef ROL_MOREAUYOSIDAPENALTY_H
#define ROL_MOREAUYOSIDAPENALTY_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::MoreauYosidaPenalty
    \brief Provides the interface to evaluate the Moreau-Yosida penalty function.

    ---
*/


namespace ROL {

template <class Real>
class MoreauYosidaPenalty : public Objective<Real> {
private:
  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<BoundConstraint<Real> > con_;

  Teuchos::RCP<Vector<Real> > l_;
  Teuchos::RCP<Vector<Real> > u_;
  Teuchos::RCP<Vector<Real> > l1_;
  Teuchos::RCP<Vector<Real> > u1_;
  Teuchos::RCP<Vector<Real> > dl1_;
  Teuchos::RCP<Vector<Real> > du1_;
  Teuchos::RCP<Vector<Real> > xlam_;
  Teuchos::RCP<Vector<Real> > v_;
  Teuchos::RCP<Vector<Real> > dv_;
  Teuchos::RCP<Vector<Real> > lam_;

  Real mu_;
  Real fval_;
  bool isConEvaluated_;
  int nfval_;
  int ngval_;

  void computePenalty(const Vector<Real> &x) {
    if ( !isConEvaluated_ ) {
      xlam_->set(x);
      xlam_->axpy(1.0/mu_,*lam_);

      l1_->set(*l_);
      l1_->axpy(-1.0,*xlam_);
      con_->pruneLowerActive(*l1_,*xlam_);
      l1_->scale(-1.0);
      l1_->plus(*l_);
      l1_->axpy(-1.0,*xlam_);

      u1_->set(*xlam_);
      u1_->axpy(-1.0,*u_);
      con_->pruneUpperActive(*u1_,*xlam_);
      u1_->scale(-1.0);
      u1_->plus(*xlam_);
      u1_->axpy(-1.0,*u_);

      dl1_->set(l1_->dual());
      con_->pruneLowerActive(*dl1_,*xlam_);
      dl1_->scale(-1.0);
      dl1_->plus(l1_->dual());

      du1_->set(u1_->dual());
      con_->pruneUpperActive(*du1_,*xlam_);
      du1_->scale(-1.0);
      du1_->plus(u1_->dual());

      isConEvaluated_ = true;
    }
  }

public:
  ~MoreauYosidaPenalty() {}

  MoreauYosidaPenalty(Objective<Real> &obj, BoundConstraint<Real> &con, 
                const ROL::Vector<Real> &x, const Real mu = 1.0)
    : mu_(mu), fval_(0.0), isConEvaluated_(false), nfval_(0), ngval_(0) {
    obj_ = Teuchos::rcp(&obj, false);
    con_ = Teuchos::rcp(&con, false);

    l_    = x.clone();
    l1_   = x.clone();
    dl1_  = x.dual().clone();
    u_    = x.clone();
    u1_   = x.clone();
    du1_  = x.dual().clone();
    xlam_ = x.clone();
    v_    = x.clone();
    dv_   = x.dual().clone();
    lam_  = x.clone();

    con_->setVectorToLowerBound(*l_);
    con_->setVectorToUpperBound(*u_);

    //lam_->zero();
    lam_->set(*u_);
    lam_->plus(*l_);
    lam_->scale(0.5);
  }

  void updateMultipliers(Real mu, const ROL::Vector<Real> &x) {
    computePenalty(x);

    lam_->set(*u1_);
    lam_->axpy(-1.0,*l1_);
    lam_->scale(mu_);

    mu_ = mu;

    nfval_ = 0;
    ngval_ = 0;
  }

  Real getObjectiveValue(void) const {
    return fval_;
  }

  int getNumberFunctionEvaluations(void) {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) {
    return ngval_;
  }

  /** \brief Update Moreau-Yosida penalty function. 

      This function updates the Moreau-Yosida penalty function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    con_->update(x,flag,iter);
    isConEvaluated_ = false;
  }

  /** \brief Compute value.

      This function returns the Moreau-Yosida penalty value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact Moreau-Yosida penalty computation.
  */
  Real value( const Vector<Real> &x, Real &tol ) {
    computePenalty(x);
    // Compute objective function value
    fval_ = obj_->value(x,tol);
    nfval_++;
    // Add value of the Moreau-Yosida penalty
    return fval_ + 0.5*mu_*(l1_->dot(*l1_) + u1_->dot(*u1_));
  }

  /** \brief Compute gradient.

      This function returns the Moreau-Yosida penalty gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact Moreau-Yosida penalty computation.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    computePenalty(x);
    // Compute gradient of objective function
    obj_->gradient(g,x,tol);
    ngval_++;
    // Add gradient of the Moreau-Yosida penalty
    g.axpy(-mu_,*dl1_);
    g.axpy(mu_,*du1_);
  }

  /** \brief Apply Hessian approximation to vector.

      This function applies the Hessian of the Moreau-Yosida penalty to the vector \f$v\f$.
      @param[out]         hv  is the the action of the Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact Moreau-Yosida penalty computation.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    computePenalty(x);
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    // Add Hessian of the Moreau-Yosida penalty
    v_->set(v);
    con_->pruneInactive(*v_,*xlam_);
    dv_->set(v_->dual());
    con_->pruneInactive(*dv_,*xlam_);
    hv.axpy(mu_,*dv_);
  }

}; // class MoreauYosidaPenalty

} // namespace ROL

#endif
