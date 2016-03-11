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

#ifndef ROL_FDIVERGENCE_HPP
#define ROL_FDIVERGENCE_HPP

#include "ROL_RiskVector.hpp"
#include "ROL_RiskMeasure.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class FDivergence : public RiskMeasure<Real> {
private:

  Real thresh_;

  Teuchos::RCP<Vector<Real> > dualVector_;

  Real xlam_;
  Real xmu_;
  Real vlam_;
  Real vmu_;

  Real valLam_;
  Real valLam2_;
  Real valMu_;
  Real valMu2_;

  bool firstReset_;

public:
  FDivergence(const Real thresh) : RiskMeasure<Real>(),
    xlam_(0), xmu_(0), vlam_(0), vmu_(0), valLam_(0), valMu_(0),
    firstReset_(true) {
    thresh_ = thresh > (Real)0 ? thresh : (Real)1.e-2;
  }

  FDivergence(Teuchos::ParameterList &parlist) : RiskMeasure<Real>(),
    xlam_(0), xmu_(0), vlam_(0), vmu_(0), valLam_(0), valMu_(0),
    firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("F-Divergence");
    Real thresh = list.get("Threshold",1.e-2);
    thresh_ = thresh > (Real)0 ? thresh : (Real)1.e-2;
  }

  virtual Real Fprimal(Real x, int deriv = 0) = 0;
  virtual Real Fdual(Real x, int deriv = 0) = 0;

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    xlam_ = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(0);
    xmu_  = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(1);
    if (firstReset_) {
      dualVector_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    dualVector_->zero();
    valLam_ = 0; valLam2_ = 0; valMu_ = 0; valMu2_ = 0;
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x, 
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0    = Teuchos::rcp_const_cast<Vector<Real> >(
            Teuchos::dyn_cast<const RiskVector<Real> >(v).getVector());
    vlam_ = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic(0);
    vmu_  = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic(1);
  }

  // Value update and get functions
  void update(const Real val, const Real weight) {
    Real r = Fdual((val-xmu_)/xlam_,0);
    RiskMeasure<Real>::val_ += weight * r;
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, gval = 0;
    sampler.sumAll(&val,&gval,1);
    return xlam_*(thresh_ + gval) + xmu_;
  }

  // Gradient update and get functions
  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real x = (val-xmu_)/xlam_;
    Real r0 = Fdual(x,0), r1 = Fdual(x,1);

    RiskMeasure<Real>::val_ += weight * r0;
    valLam_ -= weight * r1 * x;
    valMu_  -= weight * r1;

    RiskMeasure<Real>::g_->axpy(weight*r1,g);
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &gs = Teuchos::dyn_cast<RiskVector<Real> >(g);

    std::vector<Real> mygval(3), gval(3);
    mygval[0] = RiskMeasure<Real>::val_;
    mygval[1] = valLam_;
    mygval[2] = valMu_;
    sampler.sumAll(&mygval[0],&gval[0],3);

    std::vector<Real> stat(2);
    stat[0] = thresh_ + gval[0] + gval[1];
    stat[1] = (Real)1 + gval[2];
    gs.setStatistic(stat);

    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_);
    gs.setVector(*dualVector_);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv,
              const Vector<Real> &hv, const Real weight) {
    Real x = (val-xmu_)/xlam_;
    Real r1 = Fdual(x,1), r2 = Fdual(x,2);
    RiskMeasure<Real>::val_ += weight * r2 * x;
    valLam_  += weight * r2 * x * x;
    valLam2_ -= weight * r2 * gv * x;
    valMu_   += weight * r2;
    valMu2_  -= weight * r2 * gv;
    RiskMeasure<Real>::hv_->axpy(weight * r2 * (gv - vmu_ - vlam_*x)/xlam_, g);
    RiskMeasure<Real>::hv_->axpy(weight * r1, hv);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &hs = Teuchos::dyn_cast<RiskVector<Real> >(hv);

    std::vector<Real> myhval(5), hval(5);
    myhval[0] = RiskMeasure<Real>::val_;
    myhval[1] = valLam_;
    myhval[2] = valLam2_;
    myhval[3] = valMu_;
    myhval[4] = valMu2_;
    sampler.sumAll(&myhval[0],&hval[0],5);

    std::vector<Real> stat(2);
    stat[0] = (vlam_ * hval[1] + vmu_ * hval[0] + hval[2])/xlam_;
    stat[1] = (vlam_ * hval[0] + vmu_ * hval[3] + hval[4])/xlam_;
    hs.setStatistic(stat);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector_);
    hs.setVector(*dualVector_);
  }
};

}

#endif
