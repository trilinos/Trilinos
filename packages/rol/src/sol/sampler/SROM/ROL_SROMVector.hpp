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

template <class Real, class Element=Real>
class SROMVector : public Vector<Real> {
private:

  Teuchos::RCP<std::vector<Element> >  pts_vec_;
  Teuchos::RCP<std::vector<Element> >  wts_vec_;
  size_t dimension_;
  size_t numSamples_;

  Real const_pt_;
  Real const_wt_;

public:

  SROMVector(const Teuchos::RCP<std::vector<Element> > &pts_vec,
             const Teuchos::RCP<std::vector<Element> > &wts_vec,
             const Real const_pt = 1.0, const Real const_wt = 1.0)
    : pts_vec_(pts_vec), wts_vec_(wts_vec),
      const_pt_(const_pt), const_wt_(const_wt) {
    numSamples_ = wts_vec_->size();
    dimension_  = pts_vec_->size()/numSamples_; 
  }

  void set( const Vector<Real> &x ) {
    const SROMVector &ex = Teuchos::dyn_cast<const SROMVector>(x);
    for (size_t i = 0; i < numSamples_; i++) {
      for (size_t j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] = (*ex.getPoint(i))[j];
      }
      (*wts_vec_)[i] = (ex.getWeight(i)); 
    }
  }

  void plus( const Vector<Real> &x ) {
    const SROMVector &ex = Teuchos::dyn_cast<const SROMVector>(x);
    for (size_t i = 0; i < numSamples_; i++) {
      for (size_t j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] += (*ex.getPoint(i))[j];
      }
      (*wts_vec_)[i] += (ex.getWeight(i)); 
    }
  }

  void scale( const Real alpha ) {
    for (size_t i = 0; i < numSamples_; i++) {
      for (size_t j = 0; j < dimension_; j++) {
        (*pts_vec_)[i*dimension_ + j] *= alpha;
      }
      (*wts_vec_)[i] *= alpha;
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const SROMVector & ex = Teuchos::dyn_cast<const SROMVector>(x);
    Real pt_val = 0, wt_val = 0;
    for (size_t i = 0; i < numSamples_; i++) {
      for (size_t j = 0; j < dimension_; j++) {
        pt_val += (*pts_vec_)[i*dimension_ + j] * (*ex.getPoint(i))[j];
      }
      wt_val += (*wts_vec_)[i] * (ex.getWeight(i)); 
    }
    return const_pt_*pt_val + const_wt_*wt_val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    return Teuchos::rcp( new SROMVector( Teuchos::rcp(new std::vector<Element>(pts_vec_->size())),
                                         Teuchos::rcp(new std::vector<Element>(wts_vec_->size())),
                                         const_pt_, const_wt_ ) );
  }

  Teuchos::RCP<const std::vector<Element> > getPoint(const size_t i) const {
    std::vector<Element> pt(dimension_,0.);
    for (size_t j = 0; j < dimension_; j++) {
      pt[j] = (*pts_vec_)[i*dimension_ + j];
    }
    return Teuchos::rcp(new std::vector<Element>(pt));
  }

  void setPoint(const size_t i, const std::vector<Element> &pt) {
    for (size_t j = 0; j < dimension_; j++) {
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
    return dimension_;
  }

  const size_t getNumSamples(void) const {
    return numSamples_;
  }
  
}; // class SROMVector

} // namespace ROL

#endif
