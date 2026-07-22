// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEAR_OBJECTIVE_SIMOPT_H
#define ROL_LINEAR_OBJECTIVE_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"

/** @ingroup func_group
    \class ROL::LinearObjective_SimOpt
    \brief Provides the interface to evaluate linear objective functions.

    This class implements the linear objective function
    \f[
       f(x) = \langle c, x\rangle_{\mathcal{X}^*,\mathcal{X}}
    \f]
    for fixed \f$c\in\mathcal{X}^*\f$.

    ---
*/


namespace ROL {

template <class Real>
class LinearObjective_SimOpt : public Objective_SimOpt<Real> {
private:
  const Ptr<const Vector<Real>> simcost_, optcost_;

public:
  LinearObjective_SimOpt(const Ptr<const Vector<Real>> &simcost = nullPtr,
                         const Ptr<const Vector<Real>> &optcost = nullPtr)
    : simcost_(simcost), optcost_(optcost) {}

  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real valu(0), valz(0);
    if (simcost_ != nullPtr) {
      //valu = u.dot(simcost_->dual());
      valu = u.apply(*simcost_);
    }
    if (optcost_ != nullPtr) {
      //valz = z.dot(optcost_->dual());
      valz = z.apply(*optcost_);
    }
    return valu + valz;
  }

  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (simcost_ != nullPtr) {
      g.set(*simcost_);
    }
    else {
      g.zero();
    }
  }

  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (optcost_ != nullPtr) {
      g.set(*optcost_);
    }
    else {
      g.zero();
    }
  }

  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

}; // class LinearObjective_SimOpt

} // namespace ROL

#endif
