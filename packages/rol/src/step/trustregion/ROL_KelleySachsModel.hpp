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

#ifndef ROL_KELLEYSACHSMODEL_HPP
#define ROL_KELLEYSACHSMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::KelleySachsModel
    \brief Provides the interface to evaluate projected trust-region model
    functions from the Kelley-Sachs bound constrained trust-region algorithm.

    -----
*/

namespace ROL {

template<class Real>
class KelleySachsModel : public TrustRegionModel<Real> {
private:
  ROL::Ptr<BoundConstraint<Real> > bnd_;
  ROL::Ptr<Secant<Real> > secant_;
  ROL::Ptr<Vector<Real> > dual_, prim_;

  const bool useSecantPrecond_;
  const bool useSecantHessVec_;

  Real eps_;

public:

  KelleySachsModel(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                   const Vector<Real> &x, const Vector<Real> &g, const Real eps)
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      secant_(ROL::nullPtr), useSecantPrecond_(false), useSecantHessVec_(false), eps_(eps) {
    bnd_  = ROL::makePtrFromRef(bnd);
    prim_ = x.clone();
    dual_ = g.clone();
  }

  KelleySachsModel(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                   const Vector<Real> &x, const Vector<Real> &g, const Real eps,
                   const ROL::Ptr<Secant<Real> > &secant,
                   const bool useSecantPrecond, const bool useSecantHessVec)
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      secant_(secant), useSecantPrecond_(useSecantPrecond), useSecantHessVec_(useSecantHessVec), eps_(eps) {
    bnd_  = ROL::makePtrFromRef(bnd);
    prim_ = x.clone();
    dual_ = g.clone();
  }

  Real value( const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    hessVec(*dual_,s,s,tol);
    dual_->scale(static_cast<Real>(0.5));
    // Remove active components of gradient
    prim_->set(gc->dual());
    bnd_->pruneActive(*prim_,*gc,*xc,eps_);
    // Add reduced gradient to reduced hessian in direction s
    dual_->plus(prim_->dual());
    return dual_->dot(s.dual());
  }

  void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    // Apply (reduced) hessian to direction s
    hessVec(g,s,s,tol);
    // Remove active components of gradient
    prim_->set(TrustRegionModel<Real>::getGradient()->dual());
    bnd_->pruneActive(*prim_,*gc,*xc,eps_);
    // Add reduced gradient to reduced hessian in direction s
    g.plus(prim_->dual());
  }

  void hessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    // Set vnew to v
    prim_->set(v);
    // Remove elements of vnew corresponding to binding set
    bnd_->pruneActive(*prim_,*gc,*xc,eps_);
    // Apply full Hessian to reduced vector
    if ( useSecantHessVec_ ) {
      secant_->applyB(Hv,*prim_);
    }
    else {
      TrustRegionModel<Real>::getObjective()->hessVec(Hv,*prim_,*xc,tol);
    }
    // Remove elements of Hv corresponding to binding set
    bnd_->pruneActive(Hv,*gc,*xc,eps_);
    // Set vnew to v
    prim_->set(v);
    // Remove Elements of vnew corresponding to complement of binding set
    bnd_->pruneInactive(*prim_,*gc,*xc,eps_);
    dual_->set(prim_->dual());
    bnd_->pruneInactive(*dual_,*gc,*xc,eps_);
    // Fill complement of binding set elements in Hp with v
    Hv.plus(*dual_);
  }

  void invHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    // Set vnew to v
    dual_->set(v);
    // Remove elements of vnew corresponding to binding set
    bnd_->pruneActive(*dual_,*gc,*xc,eps_);
    // Apply full Hessian to reduced vector
    if ( useSecantHessVec_ ) {
      secant_->applyH(Hv,*dual_);
    }
    else {
      TrustRegionModel<Real>::getObjective()->invHessVec(Hv,*dual_,*xc,tol);
    }
    // Remove elements of Hv corresponding to binding set
    bnd_->pruneActive(Hv,*gc,*xc,eps_);
    // Set vnew to v
    dual_->set(v);
    // Remove Elements of vnew corresponding to complement of binding set
    bnd_->pruneInactive(*dual_,*gc,*xc,eps_);
    prim_->set(dual_->dual());
    bnd_->pruneInactive(*prim_,*gc,*xc,eps_);
    // Fill complement of binding set elements in Hv with v
    Hv.plus(*prim_);
  }

  void precond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    // Set vnew to v
    dual_->set(v);
    // Remove elements of vnew corresponding to binding set
    bnd_->pruneActive(*dual_,*gc,*xc,eps_);
    // Apply full Hessian to reduced vector
    if ( useSecantPrecond_ ) {
      secant_->applyH(Mv,*dual_);
    }
    else {
      TrustRegionModel<Real>::getObjective()->precond(Mv,*dual_,*xc,tol);
    }
    // Remove elements of Mv corresponding to binding set
    bnd_->pruneActive(Mv,*gc,*xc,eps_);
    // Set vnew to v
    dual_->set(v);
    // Remove Elements of vnew corresponding to complement of binding set
    bnd_->pruneInactive(*dual_,*gc,*xc,eps_);
    prim_->set(dual_->dual());
    bnd_->pruneInactive(*prim_,*gc,*xc,eps_);
    // Fill complement of binding set elements in Mv with v
    Mv.plus(*prim_);
  }

  void dualTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    // Compute T(v) = P_I(v) where P_I is the projection onto the inactive indices
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    tv.set(v);
    bnd_->pruneActive(tv,*gc,*xc,eps_);
  }

  void primalTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    // Compute T(v) = P( x + v ) - x where P is the projection onto the feasible set
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    tv.set(*xc);
    tv.plus(v);
    bnd_->project(tv);
    tv.axpy(static_cast<Real>(-1),*xc);
  }

  const ROL::Ptr<BoundConstraint<Real> > getBoundConstraint(void) const {
    return bnd_;
  }

}; 

} // namespace ROL

#endif 
