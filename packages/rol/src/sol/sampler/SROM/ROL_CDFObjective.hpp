// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CDFOBJECTIVE_H
#define ROL_CDFOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Distribution.hpp"
#include "ROL_Ptr.hpp"
#include <math.h>

namespace ROL {

template <class Real>
class CDFObjective : public Objective<Real> {
private:
  // Batch manager for parallel computation
  ROL::Ptr<BatchManager<Real> > bman_;

  // Distribution information
  std::vector<ROL::Ptr<Distribution<Real> > > dist_;
  std::vector<Real> lowerBound_;
  std::vector<Real> upperBound_;
  int dimension_;

  const Real scale_;
  const Real sqrt2_;
  const Real sqrtpi_;

  const bool optProb_;
  const bool optAtom_;

  std::vector<Real> pts_;
  std::vector<Real> wts_;

  // Number of quadrature points
  int numPoints_;

  void initializeQuadrature(void) {
    numPoints_ = 20;
    pts_.clear(); pts_.resize(numPoints_,0.);
    wts_.clear(); wts_.resize(numPoints_,0.);
    wts_[0]  = 0.1527533871307258; pts_[0]  = -0.0765265211334973;
    wts_[1]  = 0.1527533871307258; pts_[1]  =  0.0765265211334973;
    wts_[2]  = 0.1491729864726037; pts_[2]  = -0.2277858511416451;
    wts_[3]  = 0.1491729864726037; pts_[3]  =  0.2277858511416451;
    wts_[4]  = 0.1420961093183820; pts_[4]  = -0.3737060887154195;
    wts_[5]  = 0.1420961093183820; pts_[5]  =  0.3737060887154195;
    wts_[6]  = 0.1316886384491766; pts_[6]  = -0.5108670019508271;
    wts_[7]  = 0.1316886384491766; pts_[7]  =  0.5108670019508271;
    wts_[8]  = 0.1181945319615184; pts_[8]  = -0.6360536807265150;
    wts_[9]  = 0.1181945319615184; pts_[9]  =  0.6360536807265150;
    wts_[10] = 0.1019301198172404; pts_[10] = -0.7463319064601508;
    wts_[11] = 0.1019301198172404; pts_[11] =  0.7463319064601508;
    wts_[12] = 0.0832767415767048; pts_[12] = -0.8391169718222188;
    wts_[13] = 0.0832767415767048; pts_[13] =  0.8391169718222188;
    wts_[14] = 0.0626720483341091; pts_[14] = -0.9122344282513259;
    wts_[15] = 0.0626720483341091; pts_[15] =  0.9122344282513259;
    wts_[16] = 0.0406014298003869; pts_[16] = -0.9639719272779138;
    wts_[17] = 0.0406014298003869; pts_[17] =  0.9639719272779138;
    wts_[18] = 0.0176140071391521; pts_[18] = -0.9931285991850949;
    wts_[19] = 0.0176140071391521; pts_[19] =  0.9931285991850949;
    for (int i = 0; i < numPoints_; i++) {
      wts_[i] *= 0.5;
      pts_[i] += 1.; pts_[i] *= 0.5;
    }
  }

  Real valueCDF(const int dim, const Real loc,
                const ProbabilityVector<Real> &prob,
                const AtomVector<Real>        &atom) const {
    const int numSamples = prob.getNumMyAtoms();
    Real val = 0, hs = 0, xpt = 0, xwt = 0, sum = 0, half(0.5), one(1);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*atom.getAtom(k))[dim]; xwt = prob.getProbability(k);
      hs = half * (one + erf((loc-xpt)/(sqrt2_*scale_)));
      val += xwt * hs;
    }
    bman_->sumAll(&val,&sum,1);
    return sum;
  }

  Real gradientCDF(std::vector<Real> &gradx, std::vector<Real> &gradp,
             const int dim, const Real loc,
             const ProbabilityVector<Real> &prob,
             const AtomVector<Real>        &atom) const {
    const int numSamples = prob.getNumMyAtoms();
    gradx.resize(numSamples,0); gradp.resize(numSamples,0);
    Real val = 0, hs = 0, xpt = 0, xwt = 0, sum = 0, half(0.5), one(1);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*atom.getAtom(k))[dim]; xwt = prob.getProbability(k);
      hs = half * (one + erf((loc-xpt)/(sqrt2_*scale_)));
      val += xwt * hs;
      gradx[k] = -(xwt/(sqrt2_*sqrtpi_*scale_))
                 * std::exp(-std::pow((loc-xpt)/(sqrt2_*scale_),2));
      gradp[k] = hs;
    }
    bman_->sumAll(&val,&sum,1);
    return sum;
  }

  Real hessVecCDF(std::vector<Real> &hvxx, std::vector<Real> &hvxp, std::vector<Real> &hvpx,
                  std::vector<Real> &gradx, std::vector<Real> &gradp,
                  Real &sumx, Real &sump,
            const int dim, const Real loc,
             const ProbabilityVector<Real> &prob,
             const AtomVector<Real>        &atom,
             const ProbabilityVector<Real> &vprob,
             const AtomVector<Real>        &vatom) const {
    const int numSamples = prob.getNumMyAtoms();
    hvxx.resize(numSamples,0); hvxp.resize(numSamples,0); hvpx.resize(numSamples,0);
    gradx.resize(numSamples,0); gradp.resize(numSamples,0);
    sumx = 0; sump = 0;
    std::vector<Real> psum(3,0), out(3,0);
    Real val = 0, hs = 0, dval = 0, scale3 = std::pow(scale_,3);
    Real xpt = 0, xwt = 0, vpt = 0, vwt = 0, half(0.5), one(1);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*atom.getAtom(k))[dim];  xwt = prob.getProbability(k);
      vpt = (*vatom.getAtom(k))[dim]; vwt = vprob.getProbability(k);
      hs = half * (one + erf((loc-xpt)/(sqrt2_*scale_)));
      psum[0] += xwt * hs;
      dval = std::exp(-std::pow((loc-xpt)/(sqrt2_*scale_),2));
      gradx[k] = -(xwt/(sqrt2_*sqrtpi_*scale_)) * dval;
      gradp[k] = hs;
      hvxx[k] = -(xwt/(sqrt2_*sqrtpi_*scale3)) * dval * (loc-xpt) * vpt;
      hvxp[k] = -dval/(sqrt2_*sqrtpi_*scale_)*vwt;
      hvpx[k] = -dval/(sqrt2_*sqrtpi_*scale_)*vpt;
      psum[1] += vpt*gradx[k];
      psum[2] += vwt*gradp[k];
    }
    bman_->sumAll(&psum[0],&out[0],3);
    val = out[0]; sumx = out[1]; sump = out[2];
    return val;
  }

public:
  CDFObjective(const std::vector<ROL::Ptr<Distribution<Real> > > &dist,
               const ROL::Ptr<BatchManager<Real> > &bman,
               const Real scale = 1.e-2,
               const bool optProb = true, const bool optAtom = true)
    : Objective<Real>(), bman_(bman), dist_(dist), dimension_(dist.size()),
      scale_(scale), sqrt2_(std::sqrt(2.)), sqrtpi_(std::sqrt(ROL::ScalarTraits<Real>::pi())),
      optProb_(optProb), optAtom_(optAtom) {
    lowerBound_.resize(dimension_,0);
    upperBound_.resize(dimension_,0);
    for ( int i = 0; i < dimension_; i++ ) {
      lowerBound_[i] = dist[i]->lowerBound();
      upperBound_[i] = dist[i]->upperBound();
    }
    initializeQuadrature();
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const ProbabilityVector<Real> &prob = *(ex.getProbabilityVector());
    const AtomVector<Real> &atom = *(ex.getAtomVector());
    Real val(0), diff(0), pt(0), wt(0), meas(0), lb(0), one(1);
    for (int d = 0; d < dimension_; d++) {
      lb   = lowerBound_[d];
      meas = (upperBound_[d] - lb);
      meas = ((meas > ROL_EPSILON<Real>()) ? meas : one);
      for (int k = 0; k < numPoints_; k++) {
        pt = meas*pts_[k] + lb;
        wt = wts_[k]/meas;
        diff = (valueCDF(d,pt,prob,atom)-dist_[d]->evaluateCDF(pt));
        val += wt*std::pow(diff,2);
      }
    }
    return 0.5*val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.zero();
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const ProbabilityVector<Real> &prob = *(ex.getProbabilityVector());
    const AtomVector<Real> &atom = *(ex.getAtomVector());
    const int numSamples = prob.getNumMyAtoms();
    std::vector<Real> gradx(numSamples,0.), gradp(numSamples,0);
    Real diff(0), pt(0), wt(0), meas(0), lb(0), val(0), one(1);
    std::vector<Real> val_wt(numSamples,0), tmp(dimension_,0);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (int d = 0; d < dimension_; d++) {
      lb   = lowerBound_[d];
      meas = (upperBound_[d] - lb);
      meas = ((meas > ROL_EPSILON<Real>()) ? meas : one);
      for (int k = 0; k < numPoints_; k++) {
        pt = meas*pts_[k] + lb;
        wt = wts_[k]/meas;
        val = gradientCDF(gradx,gradp,d,pt,prob,atom);
        diff = (val-dist_[d]->evaluateCDF(pt));
        for (int j = 0; j < numSamples; j++) {
          (val_pt[j])[d] += wt * diff * gradx[j];
          val_wt[j]      += wt * diff * gradp[j];
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
    std::vector<Real> hvxx(numSamples,0), hvxp(numSamples,0), hvpx(numSamples,0);
    std::vector<Real> gradx(numSamples,0), gradp(numSamples,0);
    Real diff(0), pt(0), wt(0), meas(0), lb(0), val(0), sumx(0), sump(0), one(1);
    std::vector<Real> val_wt(numSamples,0), tmp(dimension_,0);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (int d = 0; d < dimension_; d++) {
      lb   = lowerBound_[d];
      meas = (upperBound_[d] - lb);
      meas = ((meas > ROL_EPSILON<Real>()) ? meas : one);
      for (int k = 0; k < numPoints_; k++) {
        pt = meas*pts_[k] + lb;
        wt = wts_[k]/meas;
        val = hessVecCDF(hvxx,hvxp,hvpx,gradx,gradp,sumx,sump,d,pt,prob,atom,vprob,vatom);
        diff = (val-dist_[d]->evaluateCDF(pt));
        for (int j = 0; j < numSamples; j++) {
          (val_pt[j])[d] += wt * ( (sump + sumx) * gradx[j] + diff * (hvxx[j] + hvxp[j]) );
          val_wt[j]      += wt * ( (sump + sumx) * gradp[j] + diff * hvpx[j] );
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
}; // class LinearCombinationObjective

} // namespace ROL

#endif
