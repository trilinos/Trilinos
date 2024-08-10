// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SHIFTEDPROXOBJECTIVE_HPP
#define ROL_SHIFTEDPROXOBJECTIVE_HPP

#include "ROL_ProxObjective.hpp"

namespace ROL {

template<typename Real>
class ShiftedProxObjective : public ProxObjective<Real> {
private:
  const Ptr<ProxObjective<Real>> prox_;
  const Ptr<Vector<Real>> x_, xtmp_;

public:
  ShiftedProxObjective(const Ptr<ProxObjective<Real>> &prox, const Ptr<Vector<Real>> &x)
    : prox_(prox), x_(x->clone()), xtmp_(x->clone()) {}

  void setX(const Vector<Real> &x) {
    x_->set(x);
  }

  Real value(const Vector<Real> x, Real &tol) override {
    xtmp_->set(*x_); xtmp_->plus(x);
    return prox_->(*xtmp_,tol);
  }

  void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) override {
    xtmp_->set(*x_); xtmp_->plus(x);
    prox_->gradient(g,*xtmp_,tol); 
  }

  void prox(Vector<Real> &x, Real gamma) override {
    x.plus(*x_);
    prox_->prox(x,gamma);
    x.axpy(static_cast<Real>(-1),*x_);
  }

};

}

#endif
