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

#ifndef ROL_RISKMEASURE_HPP
#define ROL_RISKMEASURE_HPP

#include "ROL_Vector.hpp"

namespace ROL {

template<class Real>
class RiskMeasure {
protected:
  Real val_;
  Real gv_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > hv_;

public:
  virtual ~RiskMeasure() {}

  virtual void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    this->val_ = 0.0;
    this->gv_ = 0.0;
    this->g_  = x.clone(); this->g_->zero();
    this->hv_ = x.clone(); this->hv_->zero();
    x0 = Teuchos::rcp(&const_cast<Vector<Real> &>(x),false);
  }

  virtual void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x, 
                     Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    this->val_ = 0.0;
    this->gv_ = 0.0;
    this->g_  = x.clone(); this->g_->zero();
    this->hv_ = x.clone(); this->hv_->zero();
    x0 = Teuchos::rcp(&const_cast<Vector<Real> &>(x),false);
    v0 = Teuchos::rcp(&const_cast<Vector<Real> &>(v),false);
  }

  virtual void update(const Real val, const Real weight) {
    this->val_ += weight * val;
  }

  virtual void update(const Real val, const Vector<Real> &g, const Real weight) {
    this->g_->axpy(weight,g);
  }

  virtual void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv, 
                      const Real weight) {
    this->hv_->axpy(weight,hv);
  }

  virtual Real getValue(SampleGenerator<Real> &sampler) {
    Real val = 0.0;
    sampler.sumAll(&(this->val_),&val,1);
    return val;
  }

  virtual void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*(this->g_),g);
  }

  virtual void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*(this->hv_),hv);
  }
};

}

#endif
