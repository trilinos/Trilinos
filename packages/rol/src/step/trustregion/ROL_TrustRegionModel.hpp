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

#ifndef ROL_TRUSTREGIONMODEL_H
#define ROL_TRUSTREGIONMODEL_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Secant.hpp"

/** @ingroup func_group
    \class ROL::TrustRegionModel
    \brief Provides the interface to evaluate trust-region model functions.

    ROL::TrustRegionModel provides the interface to implement a number of
    trust-region models for unconstrained and constrained optimization.
    The default implementation is the standard quadratic trust region model
    for unconstrained optimization.

    -----
*/


namespace ROL {

template <class Real>
class TrustRegionModel : public Objective<Real> {
private:
  ROL::Ptr<Objective<Real> > obj_;
  ROL::Ptr<const Vector<Real> > x_, g_;
  ROL::Ptr<Vector<Real> > dual_;
  ROL::Ptr<Secant<Real> > secant_;

  const bool useSecantPrecond_;
  const bool useSecantHessVec_;

public:

  virtual ~TrustRegionModel() {}

  TrustRegionModel(Objective<Real> &obj, const Vector<Real> &x, const Vector<Real> &g, const bool init = true)
    : secant_(ROL::nullPtr), useSecantPrecond_(false), useSecantHessVec_(false) {
    obj_  = ROL::makePtrFromRef(obj);
    x_    = ROL::makePtrFromRef(x);
    g_    = ROL::makePtrFromRef(g);
    if ( init ) {
      dual_ = g.clone();
    }
  }

  TrustRegionModel(Objective<Real> &obj, const Vector<Real> &x, const Vector<Real> &g,
                   const ROL::Ptr<Secant<Real> > &secant,
                   const bool useSecantPrecond, const bool useSecantHessVec)
    : secant_(secant), useSecantPrecond_(useSecantPrecond), useSecantHessVec_(useSecantHessVec) {
    obj_  = ROL::makePtrFromRef(obj);
    x_    = ROL::makePtrFromRef(x);
    g_    = ROL::makePtrFromRef(g);
    dual_ = g.clone();
  }

  virtual Real value( const Vector<Real> &s, Real &tol ) {
    if ( useSecantHessVec_ ) {
      secant_->applyB(*dual_,s);
    }
    else {
      obj_->hessVec(*dual_,s,*x_,tol);
    }
    dual_->scale(static_cast<Real>(0.5));
    dual_->plus(*g_);
    return dual_->dot(s.dual());
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) {
    if ( useSecantHessVec_ ) {
      secant_->applyB(g,s);
    }
    else {
      obj_->hessVec(g,s,*x_,tol);
    }
    g.plus(*g_);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    if ( useSecantHessVec_ ) {
      secant_->applyB(hv,v);
    }
    else {
      obj_->hessVec(hv,v,*x_,tol);
    }
  }

  virtual void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    if ( useSecantHessVec_ ) {
      secant_->applyH(hv,v);
    }
    else {
      obj_->invHessVec(hv,v,*x_,tol);
    }
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    if ( useSecantPrecond_ ) {
      secant_->applyH(Pv,v);
    }
    else {
      obj_->precond(Pv,v,*x_,tol);
    }
  }

  virtual void dualTransform( Vector<Real> &tv, const Vector<Real> &v ) { 
    tv.set(v);
  }

  virtual void primalTransform( Vector<Real> &tv, const Vector<Real> &v ) { 
    tv.set(v);
  }

  virtual void updatePredictedReduction(Real &pred, const Vector<Real> &s) {}

  virtual void updateActualReduction(Real &ared, const Vector<Real> &s) {}

  virtual const ROL::Ptr<const Vector<Real> > getGradient(void) const {
    return g_;
  }

  virtual const ROL::Ptr<const Vector<Real> > getIterate(void) const {
    return x_;
  }

  virtual const ROL::Ptr<Objective<Real> > getObjective(void) const {
    return obj_;
  }

  virtual const ROL::Ptr<BoundConstraint<Real> > getBoundConstraint(void) const {
    return ROL::nullPtr;
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

}; // class TrustRegionModel

} // namespace ROL


#endif
