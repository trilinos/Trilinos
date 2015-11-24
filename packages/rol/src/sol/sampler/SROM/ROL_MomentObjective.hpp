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

#ifndef ROL_MOMENTOBJECTIVE_H
#define ROL_MOMENTOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_Distribution.hpp"
#include "ROL_SROMVector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

namespace ROL {

template <class Real>
class MomentObjective : public Objective<Real> {
private:
  std::vector<std::vector<std::pair<size_t, Real> > > moments_;
  Teuchos::RCP<BatchManager<Real> > bman_;

  Real momentValue(const size_t dim, const Real power, const Real moment,
                   const PrimalSROMVector<Real> &x) const {
    const size_t numSamples = x.getNumSamples();
    Real val = 0., xpt = 0., xwt = 0., sum = 0.;
    for (size_t k = 0; k < numSamples; k++) {
      xpt = (*x.getPoint(k))[dim]; xwt = x.getWeight(k);
      val += xwt * ((power==1) ? xpt : std::pow(xpt,power));
    }
    bman_->sumAll(&val,&sum,1);
    return 0.5*std::pow((sum-moment)/moment,2);
  }

  void momentGradient(std::vector<Real> &gradx, std::vector<Real> &gradp,  Real &scale,
                const size_t dim, const Real power, const Real moment,
                const PrimalSROMVector<Real> &x) const {
    const size_t numSamples = x.getNumSamples();
    gradx.resize(numSamples,0.); gradp.resize(numSamples,0.);
    scale = 0.;
    Real xpt = 0., xwt = 0., xpow = 0., psum = 0.;
    for (size_t k = 0; k < numSamples; k++) {
      xpt = (*x.getPoint(k))[dim]; xwt = x.getWeight(k);
      xpow = ((power==1) ? 1. : ((power==2) ? xpt : std::pow(xpt,power-1)));
      psum += xwt * xpow * xpt;
      gradx[k] = xwt * xpow * power;
      gradp[k] = xpow * xpt;
    }
    bman_->sumAll(&psum,&scale,1);
    scale -= moment;
    scale /= std::pow(moment,2);
  }

  void momentHessVec(std::vector<Real> &hvx1, std::vector<Real> &hvx2, std::vector<Real> &hvx3,
                     std::vector<Real> &hvp1, std::vector<Real> &hvp2,
                     Real &scale1, Real &scale2, Real &scale3,
               const size_t dim, const Real power, const Real moment,
               const PrimalSROMVector<Real> &x, const PrimalSROMVector<Real> &v) const {
    const size_t numSamples = x.getNumSamples();
    hvx1.resize(numSamples,0.); hvx2.resize(numSamples,0.); hvx3.resize(numSamples,0.);
    hvp1.resize(numSamples,0.); hvp2.resize(numSamples,0.);
    scale1 = 0.; scale2 = 0.; scale3 = 0.;
    std::vector<Real> psum(3,0.0), scale(3,0.0);
    Real xpt = 0., xwt = 0., vpt = 0., vwt = 0.;
    Real xpow0 = 0., xpow1 = 0., xpow2 = 0.;
    const Real moment2 = std::pow(moment,2);
    for (size_t k = 0; k < numSamples; k++) {
      xpt = (*x.getPoint(k))[dim]; xwt = x.getWeight(k);
      vpt = (*v.getPoint(k))[dim]; vwt = v.getWeight(k);
      xpow2 = ((power==1) ? 0. : ((power==2) ? 1. : ((power==3) ?  xpt :
                std::pow(xpt,power-2))));
      xpow1 = ((power==1) ? 1. : xpow2 * xpt);
      xpow0 = xpow1 * xpt;
      psum[0] += xwt * xpow1 * vpt;
      psum[1] += xwt * xpow0;
      psum[2] += vwt * xpow0;
      hvx1[k] = power * xwt * xpow1;
      hvx2[k] = power * (power-1.) * xwt * xpow2 * vpt;
      hvx3[k] = power * vwt * xpow1;
      hvp1[k] = xpow0;
      hvp2[k] = power * xpow1 * vpt;
    }
    bman_->sumAll(&psum[0],&scale[0],3);
    scale1 = scale[0] * power/moment2;
    scale2 = (scale[1] - moment)/moment2 ;
    scale3 = scale[2]/moment2;
  }

public:
  MomentObjective(const std::vector<std::vector<std::pair<size_t, Real> > > &moments,
                        Teuchos::RCP<BatchManager<Real> > &bman)
    : Objective<Real>(), moments_(moments), bman_(bman) {}

  MomentObjective(const std::vector<Teuchos::RCP<Distribution<Real> > > &dist,
                  const std::vector<size_t>                             &order,
                        Teuchos::RCP<BatchManager<Real> > &bman)
    : Objective<Real>(), bman_(bman) {
    size_t numMoments = order.size();
    size_t dimension  = dist.size();
    std::vector<std::pair<size_t,Real> > data(numMoments);
    moments_.clear(); moments_.resize(dimension);
    for (size_t d = 0; d < dimension; d++) {
      for (size_t i = 0; i < numMoments; i++) {
        data[i] = std::make_pair(order[i],dist[d]->moment(order[i]));
      }
      moments_[d].assign(data.begin(),data.end());
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const PrimalSROMVector<Real> &ex = Teuchos::dyn_cast<const PrimalSROMVector<Real> >(x);
    size_t dimension  = ex.getDimension();
    Real val = 0.;
    std::vector<std::pair<size_t, Real> > data;
    for (size_t d = 0; d < dimension; d++) {
      data = moments_[d];
      for (size_t m = 0; m < data.size(); m++) {
        val += momentValue(d,(Real)data[m].first,data[m].second,ex);
      }
    }
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    DualSROMVector<Real> &eg = Teuchos::dyn_cast<DualSROMVector<Real> >(g);
    const PrimalSROMVector<Real> &ex = Teuchos::dyn_cast<const PrimalSROMVector<Real> >(x);
    size_t dimension  = ex.getDimension();
    size_t numSamples = ex.getNumSamples();
    std::vector<Real> gradx(numSamples,0.), gradp(numSamples,0.);
    Real scale = 0.;
    std::vector<std::pair<size_t, Real> > data;
    std::vector<Real> val_wt(numSamples,0.), tmp(dimension,0.);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (size_t d = 0; d < dimension; d++) {
      data = moments_[d];
      for (size_t m = 0; m < data.size(); m++) {
        momentGradient(gradx,gradp,scale,d,(Real)data[m].first,data[m].second,ex);
        for (size_t k = 0; k < numSamples; k++) {
          (val_pt[k])[d] += scale*gradx[k];
          val_wt[k]      += scale*gradp[k];
        }
      }
    }
    for (size_t k = 0; k < numSamples; k++) {
      eg.setPoint(k,val_pt[k]);
      eg.setWeight(k,val_wt[k]);
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    DualSROMVector<Real> &ehv = Teuchos::dyn_cast<DualSROMVector<Real> >(hv);
    const PrimalSROMVector<Real> &ev = Teuchos::dyn_cast<const PrimalSROMVector<Real> >(v);
    const PrimalSROMVector<Real> &ex = Teuchos::dyn_cast<const PrimalSROMVector<Real> >(x);
    const size_t dimension  = ex.getDimension();
    const size_t numSamples = ex.getNumSamples();
    std::vector<Real> hvx1(numSamples,0.), hvx2(numSamples,0.), hvx3(numSamples,0.);
    std::vector<Real> hvp1(numSamples,0.), hvp2(numSamples,0.);
    Real scale1 = 0., scale2 = 0., scale3 = 0.;
    std::vector<std::pair<size_t, Real> > data;
    std::vector<Real> val_wt(numSamples,0.), tmp(dimension,0.);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (size_t d = 0; d < dimension; d++) {
      data = moments_[d];
      for (size_t m = 0; m < data.size(); m++) {
        momentHessVec(hvx1,hvx2,hvx3,hvp1,hvp2,scale1,scale2,scale3,
                      d,(Real)data[m].first,data[m].second,ex,ev);
        for (size_t k = 0; k < numSamples; k++) {
          (val_pt[k])[d] += (scale1+scale3)*hvx1[k] + scale2*(hvx2[k]+hvx3[k]);
          val_wt[k]      += (scale1+scale3)*hvp1[k] + scale2*hvp2[k];
        }
      }
    }
    for (size_t k = 0; k < numSamples; k++) {
      ehv.setPoint(k,val_pt[k]);
      ehv.setWeight(k,val_wt[k]);
    }
  }
}; // class SROMObjective

} // namespace ROL

#endif
