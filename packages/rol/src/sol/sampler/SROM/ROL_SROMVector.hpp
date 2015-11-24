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

#ifndef ROL_SROMVECTOR_H
#define ROL_SROMVECTOR_H

#include <algorithm>
#include <cstdlib>

#include "ROL_Vector.hpp"

/** \class ROL::SROMVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real, class Element>
class SROMVector : public Vector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:

  Teuchos::RCP<BatchManager<Element> > bman_;

  Teuchos::RCP<std::vector<Element> >  pts_vec_;
  Teuchos::RCP<std::vector<Element> >  wts_vec_;

  uint dimension_;
  uint numSamples_;

  std::vector<Real> scale_wts_; // 1/pow(typw,2)
  std::vector<Real> scale_pts_; // 1/pow(typx,2)

public:

  SROMVector(const Teuchos::RCP<std::vector<Element> >  &pts_vec,
             const Teuchos::RCP<std::vector<Element> >  &wts_vec,
             const Teuchos::RCP<BatchManager<Element> > &bman,
             const std::vector<Real> &scale_pts,
             const std::vector<Real> &scale_wts)
    : bman_(bman), pts_vec_(pts_vec), wts_vec_(wts_vec) {
    numSamples_ = wts_vec_->size();
    dimension_  = 0;
    if (numSamples_ > 0) { 
      dimension_  = pts_vec_->size()/numSamples_;
    }

    Real wi = 1., xij = 1.;
    scale_wts_.clear(); scale_wts_.resize(numSamples_,1.);
    scale_pts_.clear(); scale_pts_.resize(numSamples_*dimension_,1.);
    for (uint i = 0; i < numSamples_; i++) {
      wi = std::abs(scale_wts[i]);
      if ( wi > ROL_EPSILON ) {
        scale_wts_[i] = wi;
      }
      for (uint j = 0; j < dimension_; j++) {
        xij = std::abs(scale_pts[i*dimension_ + j]);
        if ( xij > ROL_EPSILON ) {
          scale_pts_[i*dimension_ + j] = xij;
        }
      }
    }
  }

  void set( const Vector<Real> &x ) {
    const SROMVector &ex = Teuchos::dyn_cast<const SROMVector>(x);
    for (uint i = 0; i < numSamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] = (*ex.getPoint(i))[j];
      }
      (*wts_vec_)[i] = (ex.getWeight(i)); 
    }
  }

  void plus( const Vector<Real> &x ) {
    const SROMVector &ex = Teuchos::dyn_cast<const SROMVector>(x);
    for (uint i = 0; i < numSamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] += (*ex.getPoint(i))[j];
      }
      (*wts_vec_)[i] += (ex.getWeight(i)); 
    }
  }

  void scale( const Real alpha ) {
    for (uint i = 0; i < numSamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] *= alpha;
      }
      (*wts_vec_)[i] *= alpha;
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const SROMVector & ex = Teuchos::dyn_cast<const SROMVector>(x);
    Real pt_val = 0, wt_val = 0;
    for (uint i = 0; i < numSamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        pt_val += (*pts_vec_)[i*dimension_ + j]
                  * scale_pts_[i*dimension_ + j]
                  * (*ex.getPoint(i))[j];
      }
      wt_val += (*wts_vec_)[i] * scale_wts_[i] * (ex.getWeight(i)); 
    }
    // Global sum
    Real sum_pt_val = 0, sum_wt_val = 0;
    bman_->sumAll(&pt_val,&sum_pt_val,1);
    bman_->sumAll(&wt_val,&sum_wt_val,1);
    return sum_pt_val + sum_wt_val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  Teuchos::RCP<const std::vector<Element> > getPoint(const size_t i) const {
    std::vector<Element> pt(dimension_,0.);
    for (uint j = 0; j < dimension_; j++) {
      pt[j] = (*pts_vec_)[i*dimension_ + j];
    }
    return Teuchos::rcp(new std::vector<Element>(pt));
  }

  void setPoint(const size_t i, const std::vector<Element> &pt) {
    for (uint j = 0; j < dimension_; j++) {
      (*pts_vec_)[i*dimension_ + j] = pt[j];
    }
  }

  const Element getWeight(const size_t i) const {
    return (*wts_vec_)[i];
  }

  void setWeight(const size_t i, const Element wt) {
    (*wts_vec_)[i] = wt;
  }

  const size_t getDimension(void) const {
    return static_cast<size_t>(dimension_);
  }

  const size_t getNumSamples(void) const {
    return static_cast<size_t>(numSamples_);
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    for (uint i = 0; i < numSamples_; i++) {
      (*wts_vec_)[i] = f.apply((*wts_vec_)[i]);
      for (uint j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] = f.apply((*pts_vec_)[i*dimension_ + j]);
      }
    }
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const SROMVector & ex = Teuchos::dyn_cast<const SROMVector>(x);
    for (uint i = 0; i < numSamples_; i++) {
      (*wts_vec_)[i] = f.apply((*wts_vec_)[i],ex.getWeight(i));
      for (uint j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] = f.apply((*pts_vec_)[i*dimension_ + j],(*ex.getPoint(i))[j]);
      }
    }
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    for (uint i = 0; i < numSamples_; i++) {
      r.reduce((*wts_vec_)[i],result);
      for (uint j = 0; j < dimension_; j++) {
        r.reduce((*pts_vec_)[i*dimension_ + j],result);
      }
    }
    // Global sum
    Real sum = 0.;
    bman_->reduceAll(&result,&sum,r);
//    bman_->sumAll(&result,&sum,1);
    return sum;
  }

protected:

  const std::vector<Element> & getPoints(void) const {
    return *pts_vec_;
  }

  const std::vector<Element> & getWeights(void) const {
    return *wts_vec_;
  }

  const std::vector<Element> & getPointScalingVector(void) const {
    return scale_pts_;
  }

  const std::vector<Element> & getWeightScalingVector(void) const {
    return scale_wts_;
  }

  const Teuchos::RCP<BatchManager<Real> > & getBatchManager(void) const {
    return bman_;
  }
  
}; // class PrimalSROMVector

template <class Real, class Element=Real>
class PrimalSROMVector;

template <class Real, class Element=Real>
class DualSROMVector;

template <class Real, class Element>
class PrimalSROMVector : public SROMVector<Real,Element> {
  typedef typename std::vector<Real>::size_type uint;
private:

  mutable Teuchos::RCP<DualSROMVector<Real,Element> > dual_vec_;

public:

  PrimalSROMVector(const Teuchos::RCP<std::vector<Element> >  &pts_vec,
                   const Teuchos::RCP<std::vector<Element> >  &wts_vec,
                   const Teuchos::RCP<BatchManager<Element> > &bman,
                   const std::vector<Real> &scale_pts,
                   const std::vector<Real> &scale_wts)
    : SROMVector<Real,Element>(pts_vec,wts_vec,bman,scale_pts,scale_wts) {}

  Teuchos::RCP<Vector<Real> > clone() const {
    uint dimension = SROMVector<Real,Element>::getDimension();
    uint numSamples = SROMVector<Real,Element>::getNumSamples();

    return Teuchos::rcp( new PrimalSROMVector( Teuchos::rcp(new std::vector<Element>(dimension*numSamples)),
                                               Teuchos::rcp(new std::vector<Element>(numSamples)),
                                               SROMVector<Real,Element>::getBatchManager(),
                                               SROMVector<Real,Element>::getPointScalingVector(),
                                               SROMVector<Real,Element>::getWeightScalingVector() ) );
  }

  const Vector<Real> & dual() const {
    uint dimension = SROMVector<Real,Element>::getDimension();
    uint numSamples = SROMVector<Real,Element>::getNumSamples();

    std::vector<Element> tmp_pts_vec(SROMVector<Real,Element>::getPoints());
    std::vector<Element> tmp_wts_vec(SROMVector<Real,Element>::getWeights());

    std::vector<Element> scale_pts(SROMVector<Real,Element>::getPointScalingVector());
    std::vector<Element> scale_wts(SROMVector<Real,Element>::getWeightScalingVector());
    for (uint i = 0; i < numSamples; i++) {
      tmp_wts_vec[i] *= scale_wts[i];
      scale_wts[i] = 1./scale_wts[i];
      for (uint j = 0; j < dimension; j++) {
        tmp_pts_vec[i*dimension + j] *= scale_pts[i*dimension + j];
        scale_pts[i*dimension + j] = 1./scale_pts[i*dimension + j];
      }
    }
    dual_vec_ = Teuchos::rcp( new DualSROMVector<Real>( Teuchos::rcp(new std::vector<Element>(tmp_pts_vec)),
                                                        Teuchos::rcp(new std::vector<Element>(tmp_wts_vec)),
                                                        SROMVector<Real,Element>::getBatchManager(),
                                                        scale_pts, scale_wts ) );
    return *dual_vec_;
  }
  
}; // class PrimalSROMVector

template <class Real, class Element>
class DualSROMVector : public SROMVector<Real,Element> {
  typedef typename std::vector<Real>::size_type uint;
private:

  mutable Teuchos::RCP<PrimalSROMVector<Real> > dual_vec_;

public:

  DualSROMVector(const Teuchos::RCP<std::vector<Element> > &pts_vec,
                 const Teuchos::RCP<std::vector<Element> > &wts_vec,
                 const Teuchos::RCP<BatchManager<Element> > &bman,
                 const std::vector<Real> &scale_pts,
                 const std::vector<Real> &scale_wts)
    : SROMVector<Real,Element>(pts_vec,wts_vec,bman,scale_pts,scale_wts) {}

  Teuchos::RCP<Vector<Real> > clone() const {
    uint dimension = SROMVector<Real,Element>::getDimension();
    uint numSamples = SROMVector<Real,Element>::getNumSamples();

    return Teuchos::rcp( new DualSROMVector( Teuchos::rcp(new std::vector<Element>(dimension*numSamples)),
                                             Teuchos::rcp(new std::vector<Element>(numSamples)),
                                             SROMVector<Real,Element>::getBatchManager(),
                                             SROMVector<Real,Element>::getPointScalingVector(),
                                             SROMVector<Real,Element>::getWeightScalingVector() ) );
  }

  const Vector<Real> & dual() const {
    uint dimension = SROMVector<Real,Element>::getDimension();
    uint numSamples = SROMVector<Real,Element>::getNumSamples();

    std::vector<Element> tmp_pts_vec(SROMVector<Real,Element>::getPoints());
    std::vector<Element> tmp_wts_vec(SROMVector<Real,Element>::getWeights());

    std::vector<Element> scale_pts(SROMVector<Real,Element>::getPointScalingVector());
    std::vector<Element> scale_wts(SROMVector<Real,Element>::getWeightScalingVector());
    for (uint i = 0; i < numSamples; i++) {
      tmp_wts_vec[i] *= scale_wts[i];
      scale_wts[i] = 1./scale_wts[i];
      for (uint j = 0; j < dimension; j++) {
        tmp_pts_vec[i*dimension + j] *= scale_pts[i*dimension + j];
        scale_pts[i*dimension + j] = 1./scale_pts[i*dimension + j];
      }
    }
    dual_vec_ = Teuchos::rcp( new PrimalSROMVector<Real>( Teuchos::rcp(new std::vector<Element>(tmp_pts_vec)),
                                                          Teuchos::rcp(new std::vector<Element>(tmp_wts_vec)),
                                                          SROMVector<Real,Element>::getBatchManager(),
                                                          scale_pts, scale_wts ) );
    return *dual_vec_;
  }
  
  
}; // class SROMVector

} // namespace ROL

#endif
