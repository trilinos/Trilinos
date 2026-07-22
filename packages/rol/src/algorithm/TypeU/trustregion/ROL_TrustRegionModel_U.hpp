// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUSTREGIONMODEL_U_H
#define ROL_TRUSTREGIONMODEL_U_H

#include "ROL_Objective.hpp"
#include "ROL_SecantFactory.hpp"
#include "ROL_TrustRegion_U_Types.hpp"

/** @ingroup func_group
    \class ROL::TrustRegionModel_U
    \brief Provides the interface to evaluate trust-region model functions.

    ROL::TrustRegionModel_U provides the interface to implement a number of
    trust-region models for unconstrained and constrained optimization.
    The default implementation is the standard quadratic trust region model
    for unconstrained optimization.

    -----
*/

namespace ROL {

template<typename Real>
class TrustRegionModel_U : public Objective<Real> {
private:
  Ptr<Objective<Real>> obj_;
  Ptr<const Vector<Real>> x_, g_;
  Ptr<Vector<Real>> dual_;
  Real tol_;

  Ptr<Secant<Real>> secant_;
  bool useSecantPrecond_;
  bool useSecantHessVec_;

protected:
  /***************************************************************************/
  /*********  BEGIN WRAPPERS FOR HESSIAN/PRECOND APPLICATION  ****************/
  /***************************************************************************/
  void applyHessian(Vector<Real> &hv, const Vector<Real> &v, Real &tol) {
    if ( useSecantHessVec_ && secant_ != nullPtr ) {
      secant_->applyB(hv,v);
    }
    else {
      obj_->hessVec(hv,v,*x_,tol_);
    }
  }

  void applyInvHessian(Vector<Real> &hv, const Vector<Real> &v, Real &tol) {
    if ( useSecantHessVec_ && secant_ != nullPtr ) {
      secant_->applyH(hv,v);
    }
    else {
      obj_->invHessVec(hv,v,*x_,tol_);
    }
  }

  void applyPrecond(Vector<Real> &Pv, const Vector<Real> &v, Real &tol) {
    if ( useSecantPrecond_  && secant_ != nullPtr ) {
      secant_->applyH(Pv,v);
    }
    else {
      obj_->precond(Pv,v,*x_,tol_);
    }
  }
  /***************************************************************************/
  /*********  END WRAPPERS FOR HESSIAN/PRECOND APPLICATION  ******************/
  /***************************************************************************/

public:

  virtual ~TrustRegionModel_U() {}

  TrustRegionModel_U(ParameterList           &list,
                     const Ptr<Secant<Real>> &secant = nullPtr,
                     ESecantMode              mode   = SECANTMODE_BOTH)
    : obj_(nullPtr), x_(nullPtr), g_(nullPtr), secant_(secant) {
    ParameterList &slist = list.sublist("General").sublist("Secant");
    useSecantPrecond_ = slist.get("Use as Preconditioner", false);
    useSecantHessVec_ = slist.get("Use as Hessian",        false);
    if (secant_ == nullPtr) secant_ = SecantFactory<Real>(list,mode);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    dual_ = g.clone();
  }

  // Some versions of Clang will issue a warning that update hides and 
  // overloaded virtual function without this using declaration
  using Objective<Real>::update;

  void validate(Objective<Real>    &obj,
                const Vector<Real> &x,
                const Vector<Real> &g,
                ETrustRegionU       etr) {
    if ( !useSecantHessVec_ &&
        (etr == TRUSTREGION_U_DOGLEG || etr == TRUSTREGION_U_DOUBLEDOGLEG) ) {
      try {
        Real htol = std::sqrt(ROL_EPSILON<Real>());
        Ptr<Vector<Real>> v  = g.clone();
        Ptr<Vector<Real>> hv = x.clone();
        obj.invHessVec(*hv,*v,x,htol);
      }
      catch (std::exception &e) {
        useSecantHessVec_ = true;
      }
    }
  }

  virtual void setData(Objective<Real>    &obj,
                       const Vector<Real> &x,
                       const Vector<Real> &g,
                       Real &tol) {
    obj_ = makePtrFromRef(obj);
    x_   = makePtrFromRef(x);
    g_   = makePtrFromRef(g);
    tol_ = tol;
  }

  void update(const Vector<Real> &x, const Vector<Real> &s,
              const Vector<Real> &gold, const Vector<Real> &gnew,
              const Real snorm, const int iter) {
    // Update Secant Information
    if (useSecantHessVec_ || useSecantPrecond_)
      secant_->updateStorage(x,gnew,gold,s,snorm,iter);
  }

  /***************************************************************************/
  /*********  BEGIN OBJECTIVE FUNCTION DEFINITIONS  **************************/
  /***************************************************************************/
  virtual Real value( const Vector<Real> &s, Real &tol ) override {
    applyHessian(*dual_,s,tol);
    dual_->scale(static_cast<Real>(0.5));
    dual_->plus(*g_);
    return dual_->apply(s);
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) override {
    applyHessian(g,s,tol);
    g.plus(*g_);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) override {
    applyHessian(hv,v,tol);
  }

  virtual void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) override {
    applyInvHessian(hv,v,tol);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) override {
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
  /***************************************************************************/
  /*********  END ACCESSOR FUNCTIONS  ****************************************/
  /***************************************************************************/

}; // class TrustRegionModel_U

} // namespace ROL


#endif
