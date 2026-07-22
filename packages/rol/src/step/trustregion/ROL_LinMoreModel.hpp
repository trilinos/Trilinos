// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINMOREMODEL_HPP
#define ROL_LINMOREMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::LinMoreModel
    \brief Provides the interface to evaluate projected trust-region model
    functions from the Kelley-Sachs bound constrained trust-region algorithm.

    -----
*/

namespace ROL {

template<class Real>
class LinMoreModel : public TrustRegionModel<Real> {
private:
  Ptr<Vector<Real>> pwa_, dwa_;

public:

  LinMoreModel(Objective<Real> &obj, BoundConstraint<Real> &bnd,
            const Vector<Real> &x, const Vector<Real> &g,
            const Ptr<Secant<Real>> &secant = nullPtr,
            const bool useSecantPrecond = false, const bool useSecantHessVec = false)
    : TrustRegionModel<Real>::TrustRegionModel(obj,bnd,x,g,secant,useSecantPrecond,useSecantHessVec) {
    pwa_ = x.clone();
    dwa_ = g.clone();
  }

  void applyFullHessian(Vector<Real> &hv, const Vector<Real> &v, Real &tol) {
    TrustRegionModel<Real>::applyHessian(hv,v,tol);
  }

  void applyFreeHessian(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    const Real zero(0);
    pwa_->set(v);
    TrustRegionModel<Real>::getBoundConstraint()->pruneActive(*pwa_,x,zero);
    applyFullHessian(hv,*pwa_,tol);
    TrustRegionModel<Real>::getBoundConstraint()->pruneActive(hv,x,zero);
  }

  void applyFullPrecond(Vector<Real> &pv, const Vector<Real> &v, Real &tol) {
    TrustRegionModel<Real>::applyPrecond(pv,v,tol);
  }

  void applyFreePrecond(Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    const Real zero(0);
    dwa_->set(v);
    TrustRegionModel<Real>::getBoundConstraint()->pruneActive(*dwa_,x,zero);
    applyFullPrecond(pv,*dwa_,tol);
    TrustRegionModel<Real>::getBoundConstraint()->pruneActive(pv,x,zero);
  }

};

} // namespace ROL

#endif 
