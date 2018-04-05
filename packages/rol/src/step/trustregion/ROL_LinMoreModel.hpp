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
