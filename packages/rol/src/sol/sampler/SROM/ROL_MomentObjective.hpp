// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  std::vector<std::vector<std::pair<int, Real> > > moments_;
  ROL::Ptr<BatchManager<Real> > bman_;
  int dimension_;
  int numMoments_;
  const bool optProb_;
  const bool optAtom_;

  Real momentValue(const int dim, const Real power, const Real moment,
                   const ProbabilityVector<Real> &prob,
                   const AtomVector<Real>        &atom) const {
    const int numSamples = prob.getNumMyAtoms();
    Real val(0), xpt(0), xwt(0), sum(0), half(0.5), one(1), two(2);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*atom.getAtom(k))[dim]; xwt = prob.getProbability(k);
      val += xwt * ((power==one) ? xpt : std::pow(xpt,power));
    }
    bman_->sumAll(&val,&sum,1);
    Real denom = ((std::abs(moment) < ROL_EPSILON<Real>()) ? one : moment);
    return half*std::pow((sum-moment)/denom,two);
  }

  void momentGradient(std::vector<Real> &gradx, std::vector<Real> &gradp,  Real &scale,
                const int dim, const Real power, const Real moment,
                const ProbabilityVector<Real> &prob,
                const AtomVector<Real>        &atom) const {
    const int numSamples = prob.getNumMyAtoms();
    gradx.resize(numSamples,0); gradp.resize(numSamples,0);
    scale = 0;
    Real xpt(0), xwt(0), xpow(0), psum(0), one(1), two(2);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*atom.getAtom(k))[dim]; xwt = prob.getProbability(k);
      xpow = ((power==one) ? one : ((power==two) ? xpt : std::pow(xpt,power-one)));
      psum += xwt * xpow * xpt;
      gradx[k] = xwt * xpow * power;
      gradp[k] = xpow * xpt;
    }
    bman_->sumAll(&psum,&scale,1);
    scale -= moment;
    Real denom = ((std::abs(moment) < ROL_EPSILON<Real>()) ? one : moment);
    scale /= std::pow(denom,two);
  }

  void momentHessVec(std::vector<Real> &hvx1, std::vector<Real> &hvx2, std::vector<Real> &hvx3,
                     std::vector<Real> &hvp1, std::vector<Real> &hvp2,
                     Real &scale1, Real &scale2, Real &scale3,
               const int dim, const Real power, const Real moment,
               const ProbabilityVector<Real> &prob,
               const AtomVector<Real>        &atom,
               const ProbabilityVector<Real> &vprob,
               const AtomVector<Real>        &vatom) const {
    const int numSamples = prob.getNumMyAtoms();
    hvx1.resize(numSamples,0); hvx2.resize(numSamples,0); hvx3.resize(numSamples,0);
    hvp1.resize(numSamples,0); hvp2.resize(numSamples,0);
    scale1 = 0; scale2 = 0; scale3 = 0;
    std::vector<Real> psum(3,0), scale(3,0);
    Real xpt(0), xwt(0), vpt(0), vwt(0);
    Real xpow0(0), xpow1(0), xpow2(0), zero(0), one(1), two(2), three(3);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*atom.getAtom(k))[dim];  xwt = prob.getProbability(k);
      vpt = (*vatom.getAtom(k))[dim]; vwt = vprob.getProbability(k);
      xpow2 = ((power==one) ? zero : ((power==two) ? one : ((power==three) ?  xpt :
                std::pow(xpt,power-two))));
      xpow1 = ((power==one) ? one : xpow2 * xpt);
      xpow0 = xpow1 * xpt;
      psum[0] += xwt * xpow1 * vpt;
      psum[1] += xwt * xpow0;
      psum[2] += vwt * xpow0;
      hvx1[k] = power * xwt * xpow1;
      hvx2[k] = power * (power-one) * xwt * xpow2 * vpt;
      hvx3[k] = power * vwt * xpow1;
      hvp1[k] = xpow0;
      hvp2[k] = power * xpow1 * vpt;
    }
    bman_->sumAll(&psum[0],&scale[0],3);
    Real denom = ((std::abs(moment) < ROL_EPSILON<Real>()) ? one : moment);
    Real denom2 = denom*denom;
    //const Real moment2 = std::pow(moment,2);
    scale1 = scale[0] * power/denom2;
    scale2 = (scale[1] - moment)/denom2 ;
    scale3 = scale[2]/denom2;
  }

public:
  MomentObjective(const std::vector<std::vector<std::pair<int, Real> > > &moments,
                  const ROL::Ptr<BatchManager<Real> > &bman,
                  const bool optProb = true, const bool optAtom = true)
    : Objective<Real>(), moments_(moments), bman_(bman),
      optProb_(optProb), optAtom_(optAtom) {
    dimension_ = moments_.size();
    numMoments_ = moments_[0].size();
  }

  MomentObjective(const std::vector<ROL::Ptr<Distribution<Real> > > &dist,
                  const std::vector<int>                             &order,
                  const ROL::Ptr<BatchManager<Real> > &bman,
                  const bool optProb = true, const bool optAtom = true)
    : Objective<Real>(), bman_(bman), optProb_(optProb), optAtom_(optAtom) {
    numMoments_ = order.size();
    dimension_  = dist.size();
    std::vector<std::pair<int,Real> > data(numMoments_);
    moments_.clear(); moments_.resize(dimension_);
    for (int d = 0; d < dimension_; d++) {
      for (int i = 0; i < numMoments_; i++) {
        data[i] = std::make_pair(order[i],dist[d]->moment(order[i]));
      }
      moments_[d].assign(data.begin(),data.end());
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const ProbabilityVector<Real> &prob = *(ex.getProbabilityVector());
    const AtomVector<Real> &atom = *(ex.getAtomVector());
    Real val(0);
    std::vector<std::pair<int, Real> > data;
    for (int d = 0; d < dimension_; d++) {
      data = moments_[d];
      for (int m = 0; m < numMoments_; m++) {
        val += momentValue(d,(Real)data[m].first,data[m].second,prob,atom);
      }
    }
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.zero();
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const ProbabilityVector<Real> &prob = *(ex.getProbabilityVector());
    const AtomVector<Real> &atom = *(ex.getAtomVector());
    int numSamples = prob.getNumMyAtoms();
    std::vector<Real> gradx(numSamples,0), gradp(numSamples,0);
    Real scale(0);
    std::vector<std::pair<int, Real> > data;
    std::vector<Real> val_wt(numSamples,0), tmp(dimension_,0);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (int d = 0; d < dimension_; d++) {
      data = moments_[d];
      for (int m = 0; m < numMoments_; m++) {
        momentGradient(gradx,gradp,scale,d,(Real)data[m].first,data[m].second,prob,atom);
        for (int k = 0; k < numSamples; k++) {
          (val_pt[k])[d] += scale*gradx[k];
          val_wt[k]      += scale*gradp[k];
        }
      }
    }
    SROMVector<Real> &eg = dynamic_cast<SROMVector<Real>&>(g);
    ProbabilityVector<Real> &gprob = *(eg.getProbabilityVector());
    AtomVector<Real> &gatom = *(eg.getAtomVector());
    for (int k = 0; k < numSamples; k++) {
      if ( optProb_ ) {
        gprob.setProbability(k,val_wt[k]);
      }
      if ( optAtom_ ) {
        gatom.setAtom(k,val_pt[k]);
      }
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.zero();
    const SROMVector<Real> &ev = dynamic_cast<const SROMVector<Real>&>(v);
    const ProbabilityVector<Real> &vprob = *(ev.getProbabilityVector());
    const AtomVector<Real> &vatom = *(ev.getAtomVector());
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const ProbabilityVector<Real> &prob = *(ex.getProbabilityVector());
    const AtomVector<Real> &atom = *(ex.getAtomVector());
    const int numSamples = prob.getNumMyAtoms();
    std::vector<Real> hvx1(numSamples,0), hvx2(numSamples,0), hvx3(numSamples,0);
    std::vector<Real> hvp1(numSamples,0), hvp2(numSamples,0);
    Real scale1(0), scale2(0), scale3(0);
    std::vector<std::pair<int, Real> > data;
    std::vector<Real> val_wt(numSamples,0), tmp(dimension_,0);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (int d = 0; d < dimension_; d++) {
      data = moments_[d];
      for (int m = 0; m < numMoments_; m++) {
        momentHessVec(hvx1,hvx2,hvx3,hvp1,hvp2,scale1,scale2,scale3,
                      d,(Real)data[m].first,data[m].second,prob,atom,vprob,vatom);
        for (int k = 0; k < numSamples; k++) {
          (val_pt[k])[d] += (scale1+scale3)*hvx1[k] + scale2*(hvx2[k]+hvx3[k]);
          val_wt[k]      += (scale1+scale3)*hvp1[k] + scale2*hvp2[k];
        }
      }
    }
    SROMVector<Real> &ehv = dynamic_cast<SROMVector<Real>&>(hv);
    ProbabilityVector<Real> &hprob = *(ehv.getProbabilityVector());
    AtomVector<Real> &hatom = *(ehv.getAtomVector());
    for (int k = 0; k < numSamples; k++) {
      if ( optProb_ ) {
        hprob.setProbability(k,val_wt[k]);
      }
      if ( optAtom_ ) {
        hatom.setAtom(k,val_pt[k]);
      }
    }
  }
}; // class SROMObjective

} // namespace ROL

#endif
