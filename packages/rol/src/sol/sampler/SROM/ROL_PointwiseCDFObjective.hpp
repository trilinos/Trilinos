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

#ifndef ROL_POINTWISECDFOBJECTIVE_H
#define ROL_POINTWISECDFOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Distribution.hpp"
#include "ROL_Ptr.hpp"
#include <math.h>

namespace ROL {

template <class Real>
class PointwiseCDFObjective : public Objective<Real> {
private:
  std::vector<ROL::Ptr<Distribution<Real> > > dist_;
  ROL::Ptr<BatchManager<Real> > bman_;
  const Real scale_;
  const Real sqrt2_;
  const Real sqrtpi_;

  Real valueCDF(const int dim, const Real loc, const SROMVector<Real> &x) const {
    const int numSamples = x.getNumSamples();
    Real val = 0., hs = 0., xpt = 0., xwt = 0.;
    for (int k = 0; k < numSamples; k++) {
      xpt = (*x.getPoint(k))[dim]; xwt = x.getWeight(k);
      hs = 0.5 * (1. + erf((loc-xpt)/(sqrt2_*scale_)));
      val += xwt * hs;
    }
    return val;
  }

  Real gradientCDF(std::vector<Real> &gradx, std::vector<Real> &gradp,
             const int dim, const Real loc, const SROMVector<Real> &x) const {
    const int numSamples = x.getNumSamples();
    gradx.resize(numSamples,0.); gradp.resize(numSamples,0.);
    Real val = 0., hs = 0., xpt = 0., xwt = 0.;
    for (int k = 0; k < numSamples; k++) {
      xpt = (*x.getPoint(k))[dim]; xwt = x.getWeight(k);
      hs = 0.5 * (1. + erf((loc-xpt)/(sqrt2_*scale_)));
      val += xwt * hs;
      gradx[k] = -(xwt/(sqrt2_*sqrtpi_*scale_))
                 * std::exp(-std::pow((loc-xpt)/(sqrt2_*scale_),2));
      gradp[k] = hs;
    }
    return val;
  }

  Real hessVecCDF(std::vector<Real> &hvx,
            const int dim, const Real loc, const SROMVector<Real> &x, const SROMVector<Real> &v) const {
    const int numSamples = x.getNumSamples();
    hvx.resize(numSamples,0.);
    Real val = 0., hs = 0., xpt = 0., xwt = 0., scale3 = std::pow(scale_,3);
    for (int k = 0; k < numSamples; k++) {
      xpt = (*x.getPoint(k))[dim]; xwt = x.getWeight(k);
      hs = 0.5 * (1. + erf((loc-xpt)/(sqrt2_*scale_)));
      val += xwt * hs;
      hvx[k] = -(xwt/(sqrt2_*sqrtpi_*scale3))
               * std::exp(-std::pow((loc-xpt)/(sqrt2_*scale_),2)) * (loc-xpt);
    }
    return val;
  }

public:
  PointwiseCDFObjective(const std::vector<ROL::Ptr<Distribution<Real> > > &dist,
                              ROL::Ptr<BatchManager<Real> >               &bman,
                        const Real scale = 1.e-2)
    : Objective<Real>(), dist_(dist), bman_(bman), scale_(scale),
      sqrt2_(std::sqrt(2.)), sqrtpi_(std::sqrt(Teuchos::ScalarTraits<Real>::pi())) {}

  Real value( const Vector<Real> &x, Real &tol ) {
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const int dimension  = ex.getDimension();
    const int numSamples = ex.getNumSamples();
    Real val = 0., diff = 0., xpt = 0., sum = 0.;
    for (int d = 0; d < dimension; d++) {
      for (int k = 0; k < numSamples; k++) {
        xpt = (*ex.getPoint(k))[d];
        diff = (valueCDF(d,xpt,ex)-dist_[d]->evaluateCDF(xpt));
        val += std::pow(diff,2);
      }
    }
    bman_->sumAll(&val,&sum,1);
    return 0.5*sum;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    SROMVector<Real> &eg = dynamic_cast<SROMVector<Real>&>(g);
    const SROMVector<Real> &ex = dynamic_cast<const SROMVector<Real>&>(x);
    const int dimension  = ex.getDimension();
    const int numSamples = ex.getNumSamples();
    std::vector<Real> gradx(numSamples,0.), gradp(numSamples,0.);
    Real diff = 0., xpt = 0., val = 0., sum = 0.;
    std::vector<Real> val_wt(numSamples,0.), tmp(dimension,0.);
    std::vector<std::vector<Real> > val_pt(numSamples,tmp);
    for (int d = 0; d < dimension; d++) {
      for (int k = 0; k < numSamples; k++) {
        xpt = (*ex.getPoint(k))[d];
        val = gradientCDF(gradx,gradp,d,xpt,ex);
        diff = (val-dist_[d]->evaluateCDF(xpt));
        sum = 0.;
        for (int j = 0; j < numSamples; j++) {
          (val_pt[j])[d] += diff * gradx[j];
          val_wt[j]      += diff * gradp[j];
          sum            -= gradx[j]; 
        }
        (val_pt[k])[d] += diff * (sum - dist_[d]->evaluatePDF(xpt));
      }
    }
    for (int k = 0; k < numSamples; k++) {
      eg.setPoint(k,val_pt[k]);
      eg.setWeight(k,val_wt[k]);
    }
  }

//  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
//  }
}; // class LinearCombinationObjective

} // namespace ROL

#endif
