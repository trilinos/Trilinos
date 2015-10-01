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

/** \file
    \brief  Contains definitions for std::vector bound constraints.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_SROMBOUNDCONSTRAINT_HPP
#define ROL_SROMBOUNDCONSTRAINT_HPP

#include "ROL_SROMVector.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {

template<class Real>
class SROMBoundConstraint : public BoundConstraint<Real> {
private:
  std::vector<Real> pt_lo_;
  std::vector<Real> pt_up_;

  size_t dimension_;

  Real min_diff_;
  Real scale_;

public:
  
  SROMBoundConstraint(const std::vector<Real> &pt_lo,
                      const std::vector<Real> &pt_up,
                      const Real scale = 1.0)
    : pt_lo_(pt_lo), pt_up_(pt_up), scale_(scale) {
    dimension_ = pt_lo_.size();
    min_diff_  = 1.0;
    Real diff  = 0.0;
    for ( size_t i = 0; i < dimension_; i++ ) {
      diff = pt_up[i] - pt_lo_[i];
      min_diff_ = ( (min_diff_ < diff) ? min_diff_ : diff );
    }
    min_diff_ *= 0.5;
  }

  SROMBoundConstraint(const size_t dimension, const Real scale = 1.0)
    : dimension_(dimension), min_diff_(1.0), scale_(scale) {
    pt_lo_.clear(); pt_lo_.assign(dimension_,-0.1*ROL_OVERFLOW);
    pt_up_.clear(); pt_up_.assign(dimension_, 0.1*ROL_OVERFLOW);
  }

  bool isFeasible( const Vector<Real> &x ) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    size_t cnt = 1;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      for (size_t j = 0; j < dimension_; j++) {
        cnt *= ( ( ((*ex.getPoint(i))[j] >= pt_lo_[j]) && ((*ex.getPoint(i))[j] <= pt_up_[j]) ) ? 1 : 0 );
      }
      cnt *= ( ((ex.getWeight(i) >= 0.0) && (ex.getWeight(i) <= 1.0)) ? 1 : 0 );
    }
    return ((cnt==0) ? false : true);
  }

  void project( Vector<Real> &x ) {
    SROMVector<Real> &ex = Teuchos::dyn_cast<SROMVector<Real> >(x);
    std::vector<Real> pt(dimension_,0.0);
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = std::max(pt_lo_[j],std::min(pt_up_[j],(*ex.getPoint(i))[j]));
      }
      ex.setPoint(i,pt);
      ex.setWeight(i,std::max(0.0,std::min(1.0,ex.getWeight(i))));
    }
  }

  void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    SROMVector<Real> &ev = Teuchos::dyn_cast<SROMVector<Real> >(v);
    Real epsn = std::min(scale_*eps,min_diff_);
    std::vector<Real> pt(dimension_,0.0);
    Teuchos::RCP<const std::vector<Real> > xpt, vpt;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      xpt = ex.getPoint(i);
      vpt = ev.getPoint(i);
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = ( ((*xpt)[j] <= pt_lo_[j]+epsn) ? 0.0 : (*vpt)[j] );
      }
      ev.setPoint(i,pt);
      if ( ex.getWeight(i) <= 0.0+epsn ) {
        ev.setWeight(i,0.0);
      }
    }
  }

  void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    SROMVector<Real> &ev = Teuchos::dyn_cast<SROMVector<Real> >(v);
    Real epsn = std::min(scale_*eps,min_diff_);
    std::vector<Real> pt(dimension_,0.0);
    Teuchos::RCP<const std::vector<Real> > xpt, vpt;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      xpt = ex.getPoint(i);
      vpt = ev.getPoint(i);
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = ( ((*xpt)[j] >= pt_up_[j]-epsn) ? 0.0 : (*vpt)[j] );
      }
      ev.setPoint(i,pt);
      if ( ex.getWeight(i) >= 1.0-epsn ) {
        ev.setWeight(i,0.0);
      }
    }
  }

  void pruneActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    SROMVector<Real> &ev = Teuchos::dyn_cast<SROMVector<Real> >(v);
    Real epsn = std::min(scale_*eps,min_diff_);
    std::vector<Real> pt(dimension_,0.0);
    Teuchos::RCP<const std::vector<Real> > xpt, vpt;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      xpt = ex.getPoint(i);
      vpt = ev.getPoint(i);
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = ( ((*xpt)[j] >= pt_up_[j]-epsn || (*xpt)[j] <= pt_lo_[j]+epsn) ? 0.0 : (*vpt)[j] );
      }
      ev.setPoint(i,pt);
      if ( ex.getWeight(i) >= 1.0-epsn || ex.getWeight(i) <= 0.0+epsn ) {
        ev.setWeight(i,0.0);
      }
    }
  }

  void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    const SROMVector<Real> &eg = Teuchos::dyn_cast<const SROMVector<Real> >(g);
    SROMVector<Real> &ev = Teuchos::dyn_cast<SROMVector<Real> >(v);
    Real epsn = std::min(scale_*eps,min_diff_);
    std::vector<Real> pt(dimension_,0.0);
    Teuchos::RCP<const std::vector<Real> > xpt, vpt, gpt;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      xpt = ex.getPoint(i);
      gpt = eg.getPoint(i);
      vpt = ev.getPoint(i);
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = ( ((*xpt)[j] <= pt_lo_[j]+epsn && (*gpt)[i] > 0.0) ? 0.0 : (*vpt)[j] );
      }
      ev.setPoint(i,pt);
      if ( ex.getWeight(i) <= 0.0+epsn && eg.getWeight(i) > 0.0 ) {
        ev.setWeight(i,0.0);
      }
    }
  }

  void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    const SROMVector<Real> &eg = Teuchos::dyn_cast<const SROMVector<Real> >(g);
    SROMVector<Real> &ev = Teuchos::dyn_cast<SROMVector<Real> >(v);
    Real epsn = std::min(scale_*eps,min_diff_);
    std::vector<Real> pt(dimension_,0.0);
    Teuchos::RCP<const std::vector<Real> > xpt, vpt, gpt;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      xpt = ex.getPoint(i);
      gpt = eg.getPoint(i);
      vpt = ev.getPoint(i);
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = ( ((*xpt)[j] >= pt_up_[j]-epsn && (*gpt)[i] < 0.0) ? 0.0 : (*vpt)[j] );
      }
      ev.setPoint(i,pt);
      if ( ex.getWeight(i) >= 1.0-epsn && eg.getWeight(i) < 0.0 ) {
        ev.setWeight(i,0.0);
      }
    }
  }

  void pruneActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(x);
    const SROMVector<Real> &eg = Teuchos::dyn_cast<const SROMVector<Real> >(g);
    SROMVector<Real> &ev = Teuchos::dyn_cast<SROMVector<Real> >(v);
    Real epsn = std::min(scale_*eps,min_diff_);
    std::vector<Real> pt(dimension_,0.0);
    Teuchos::RCP<const std::vector<Real> > xpt, vpt, gpt;
    for (size_t i = 0; i < ex.getNumSamples(); i++) {
      xpt = ex.getPoint(i);
      gpt = eg.getPoint(i);
      vpt = ev.getPoint(i);
      for (size_t j = 0; j < dimension_; j++) {
        pt[j] = ( (((*xpt)[j] >= pt_up_[j]-epsn && (*gpt)[i] < 0.0) ||
                   ((*xpt)[j] <= pt_lo_[j]+epsn && (*gpt)[i] > 0.0)) ? 0.0 : (*vpt)[j] );
      }
      ev.setPoint(i,pt);
      if ( (ex.getWeight(i) >= 1.0-epsn && eg.getWeight(i) < 0.0) ||
           (ex.getWeight(i) <= 0.0+epsn && eg.getWeight(i) > 0.0) ) {
        ev.setWeight(i,0.0);
      }
    }
  }

  void setVectorToUpperBound( ROL::Vector<Real> &u ) {
    SROMVector<Real> &eu = Teuchos::dyn_cast<SROMVector<Real> >(u);
    for (size_t i = 0; i < eu.getNumSamples(); i++) {
      eu.setPoint(i,pt_up_);
      eu.setWeight(i,1.0);
    }
  }

  void setVectorToLowerBound( ROL::Vector<Real> &l ) {
    SROMVector<Real> &el = Teuchos::dyn_cast<SROMVector<Real> >(l);
    for (size_t i = 0; i < el.getNumSamples(); i++) {
      el.setPoint(i,pt_lo_);
      el.setWeight(i,0.0);
    }
  }
};

}// End ROL Namespace

#endif
