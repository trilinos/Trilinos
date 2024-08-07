// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
