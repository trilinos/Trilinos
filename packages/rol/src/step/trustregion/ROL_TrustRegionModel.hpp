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
  Ptr<Objective<Real>> obj_;
  Ptr<BoundConstraint<Real>> bnd_;
  Ptr<const Vector<Real>> x_, g_;
  Ptr<Vector<Real>> dual_;
  Ptr<Secant<Real>> secant_;

  const bool useSecantPrecond_;
  const bool useSecantHessVec_;

  bool init_;

  void initialize(const Vector<Real> &s) {
    if (!init_) {
      dual_ = s.dual().clone();
      init_ = true;
    }
  }

protected:
  /***************************************************************************/
  /*********  BEGIN WRAPPERS FOR HESSIAN/PRECOND APPLICATION  ****************/
  /***************************************************************************/
  void applyHessian(Vector<Real> &hv, const Vector<Real> &v, Real &tol) {
    if ( useSecantHessVec_ && secant_ != nullPtr ) {
      secant_->applyB(hv,v);
    }
    else {
      obj_->hessVec(hv,v,*x_,tol);
    }
  }

  void applyInvHessian(Vector<Real> &hv, const Vector<Real> &v, Real &tol) {
    if ( useSecantHessVec_ && secant_ != nullPtr ) {
      secant_->applyH(hv,v);
    }
    else {
      obj_->invHessVec(hv,v,*x_,tol);
    }
  }

  void applyPrecond(Vector<Real> &Pv, const Vector<Real> &v, Real &tol) {
    if ( useSecantPrecond_  && secant_ != nullPtr ) {
      secant_->applyH(Pv,v);
    }
    else {
      obj_->precond(Pv,v,*x_,tol);
    }
  }
  /***************************************************************************/
  /*********  END WRAPPERS FOR HESSIAN/PRECOND APPLICATION  ******************/
  /***************************************************************************/

public:

  virtual ~TrustRegionModel() {}

  TrustRegionModel(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                   const Vector<Real> &x, const Vector<Real> &g,
                   const Ptr<Secant<Real>> &secant = nullPtr,
                   const bool useSecantPrecond = false, const bool useSecantHessVec = false)
    : obj_(makePtrFromRef(obj)), bnd_(makePtrFromRef(bnd)),
      x_(makePtrFromRef(x)), g_(makePtrFromRef(g)),
      secant_(secant), useSecantPrecond_(useSecantPrecond), useSecantHessVec_(useSecantHessVec),
      init_(false) {}

  // Some versions of Clang will issue a warning that update hides and 
  // overloaded virtual function without this using declaration
  using Objective<Real>::update;

  virtual void update(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                const Vector<Real> &x, const Vector<Real> &g,
                const Ptr<Secant<Real>> &secant = nullPtr) {
    obj_    = makePtrFromRef(obj);
    bnd_    = makePtrFromRef(bnd);
    x_      = makePtrFromRef(x);
    g_      = makePtrFromRef(g);
    secant_ = secant;
  }

  /***************************************************************************/
  /*********  BEGIN OBJECTIVE FUNCTION DEFINITIONS  **************************/
  /***************************************************************************/
  virtual Real value( const Vector<Real> &s, Real &tol ) {
    initialize(s);
    applyHessian(*dual_,s,tol);
    dual_->scale(static_cast<Real>(0.5));
    dual_->plus(*g_);
    return dual_->dot(s.dual());
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) {
    applyHessian(g,s,tol);
    g.plus(*g_);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    applyHessian(hv,v,tol);
  }

  virtual void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    applyInvHessian(hv,v,tol);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    applyPrecond(Pv,v,tol);
  }
  /***************************************************************************/
  /*********  END OBJECTIVE FUNCTION DEFINITIONS  ****************************/
  /***************************************************************************/

  /***************************************************************************/
  /*********  BEGIN ACCESSOR FUNCTIONS  **************************************/
  /***************************************************************************/
  virtual const Ptr<const Vector<Real>> getGradient(void) const {
    return g_;
  }

  virtual const Ptr<const Vector<Real>> getIterate(void) const {
    return x_;
  }

  virtual const Ptr<Objective<Real>> getObjective(void) const {
    return obj_;
  }

  virtual const Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    if (!bnd_->isActivated()) {
      return nullPtr;
    }
    return bnd_;
  }
  /***************************************************************************/
  /*********  END ACCESSOR FUNCTIONS  ****************************************/
  /***************************************************************************/

  virtual void dualTransform( Vector<Real> &tv, const Vector<Real> &v ) { 
    tv.set(v);
  }

  virtual void primalTransform( Vector<Real> &tv, const Vector<Real> &v ) { 
    tv.set(v);
  }

  virtual void updatePredictedReduction(Real &pred, const Vector<Real> &s) {}

  virtual void updateActualReduction(Real &ared, const Vector<Real> &s) {}

}; // class TrustRegionModel

} // namespace ROL


#endif
