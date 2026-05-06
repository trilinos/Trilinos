// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_AUGMENTEDOBJECTIVE_HPP
#define ROL_OED_AUGMENTEDOBJECTIVE_HPP

#include <vector>
#include "ROL_Ptr.hpp"
#include "ROL_Objective.hpp"
#include "ROL_OED_DesignVector.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class AugmentedObjective : public Objective<Real> {
private:
  const std::vector<Ptr<Objective<Real>>> penVec_;
  const Ptr<Objective<Real>> robObj_;
  const std::vector<Real> weights_;
  Ptr<Vector<Real>> dwa_;

  void initialize(const DesignVector<Real> &d) {
  }

public:
  AugmentedObjective(const std::vector<Ptr<Objective<Real>>> &penVec,
                     const std::vector<Real> &weights)
    : penVec_(penVec), robObj_(makePtr<RobustObjective<Real>>()), weights_(weights),
      dwa_(nullPtr) {}

  void update(const Vector<Real> &x, UpdateType type, int iter = -1 ) override {
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    robObj_->update(*xs.getVector(),type,iter);
    for (unsigned i = 0u; i < weights_.size(); ++i)
      penVec_[i]->update(*xs.getVector(),type,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) override {
    const DesignVector<Real>& xs = static_cast<const DesignVector<Real>&>(x);
    Real val = robObj_->value(x,tol);
    for (unsigned i = 0u; i < weights_.size(); ++i)
      val += weights_[i] * penVec_[i]->value(*xs.getVector(),tol);
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override {
    DesignVector<Real>& gs = static_cast<DesignVector<Real>&>(g);
    const DesignVector<Real>& xs = static_cast<const DesignVector<Real>&>(x);
    if (dwa_ == nullPtr) dwa_ = gs.getVector()->clone();
    robObj_->gradient(g,x,tol);
    for (unsigned i = 0u; i < weights_.size(); ++i) {
      penVec_[i]->gradient(*dwa_,*xs.getVector(),tol);
      gs.getVector()->axpy(weights_[i],*dwa_);
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override {
    DesignVector<Real>& hs = static_cast<DesignVector<Real>&>(hv);
    const DesignVector<Real>& vs = static_cast<const DesignVector<Real>&>(v);
    const DesignVector<Real>& xs = static_cast<const DesignVector<Real>&>(x);
    if (dwa_ == nullPtr) dwa_ = hs.getVector()->clone();
    robObj_->hessVec(hv,v,x,tol);
    for (unsigned i = 0u; i < weights_.size(); ++i) {
      penVec_[i]->hessVec(*dwa_,*vs.getVector(),*xs.getVector(),tol);
      hs.getVector()->axpy(weights_[i],*dwa_);
    }
  }

}; // class AugmentedObjective

} // End OED Namespace
} // End ROL Namespace

#endif
