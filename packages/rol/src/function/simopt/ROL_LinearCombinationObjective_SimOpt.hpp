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

#ifndef ROL_LINEARCOMBINATIONOBJECTIVE_SIMOPT_H
#define ROL_LINEARCOMBINATIONOBJECTIVE_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {

template <class Real>
class LinearCombinationObjective_SimOpt : public Objective_SimOpt<Real> {
private:
  const std::vector<ROL::Ptr<Objective_SimOpt<Real> > > obj_;
  std::vector<Real> weights_;
  size_t size_;

  ROL::Ptr<Vector<Real> > udual_, zdual_;
  bool uinitialized_, zinitialized_;

public:
  LinearCombinationObjective_SimOpt(const std::vector<ROL::Ptr<Objective_SimOpt<Real> > > &obj)
    : Objective_SimOpt<Real>(), obj_(obj),
      udual_(ROL::nullPtr), zdual_(ROL::nullPtr),
      uinitialized_(false), zinitialized_(false) {
    size_ = obj_.size();
    weights_.clear(); weights_.assign(size_,static_cast<Real>(1));
  }

  LinearCombinationObjective_SimOpt(const std::vector<Real> &weights,
                                    const std::vector<ROL::Ptr<Objective_SimOpt<Real> > > &obj)
    : Objective_SimOpt<Real>(), obj_(obj),
      weights_(weights), size_(weights.size()),
      udual_(ROL::nullPtr), zdual_(ROL::nullPtr),
      uinitialized_(false), zinitialized_(false) {}

  void update(const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1) {
    for (size_t i=0; i<size_; ++i) {
      obj_[i]->update(u,z,flag,iter);
    }
  }

  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real val(0);
    for (size_t i = 0; i < size_; ++i) {
      val += weights_[i]*obj_[i]->value(u,z,tol);
    } 
    return val;
  }

  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (!uinitialized_) {
      udual_ = g.clone();
      uinitialized_ = true;
    }
    g.zero();
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->gradient_1(*udual_,u,z,tol);
      g.axpy(weights_[i],*udual_);
    }
  }

  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (!zinitialized_) {
      zdual_ = g.clone();
      zinitialized_ = true;
    }
    g.zero();
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->gradient_2(*zdual_,u,z,tol);
      g.axpy(weights_[i],*zdual_);
    }
  }

  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (!uinitialized_) {
      udual_ = hv.clone();
      uinitialized_ = true;
    }
    hv.zero();
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->hessVec_11(*udual_,v,u,z,tol);
      hv.axpy(weights_[i],*udual_);
    }
  }

  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (!uinitialized_) {
      udual_ = hv.clone();
      uinitialized_ = true;
    }
    hv.zero();
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->hessVec_12(*udual_,v,u,z,tol);
      hv.axpy(weights_[i],*udual_);
    }
  }

  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (!zinitialized_) {
      zdual_ = hv.clone();
      zinitialized_ = true;
    }
    hv.zero();
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->hessVec_21(*zdual_,v,u,z,tol);
      hv.axpy(weights_[i],*zdual_);
    }
  }

  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (!zinitialized_) {
      zdual_ = hv.clone();
      zinitialized_ = true;
    }
    hv.zero();
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->hessVec_22(*zdual_,v,u,z,tol);
      hv.axpy(weights_[i],*zdual_);
    }
  }

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) {
    Objective_SimOpt<Real>::setParameter(param);
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->setParameter(param);
    }
  }
}; // class LinearCombinationObjective

} // namespace ROL

#endif
