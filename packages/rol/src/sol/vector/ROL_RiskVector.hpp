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

#ifndef ROL_RISKVECTOR_HPP
#define ROL_RISKVECTOR_HPP

#include "ROL_Vector.hpp"
#include "ROL_RiskMeasureInfo.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace ROL {

template<class Real> 
class RiskVector : public Vector<Real> {
  typedef typename std::vector<Real>::size_type uint;

private:
  std::vector<Real> stat_;
  Teuchos::RCP<Vector<Real> > vec_;
  bool augmented_;
  uint nStat_;

  mutable Teuchos::RCP<Vector<Real> > dual_vec1_;
  mutable Teuchos::RCP<RiskVector<Real> > dual_vec_;

public:
  RiskVector( Teuchos::ParameterList &parlist,
        const Teuchos::RCP<Vector<Real> > &vec,
        const Real stat = 1 )
    : vec_(vec), augmented_(false), nStat_(0) {
    // Initialize dual vector storage
    dual_vec1_ = vec->dual().clone();
    // Get risk measure information
    std::string name;
    std::vector<Real> lower, upper;
    bool activated(false);
    int nStat(0);
    RiskMeasureInfo<Real>(parlist,name,nStat,lower,upper,activated);
    augmented_ = (nStat > 0) ? true : false;
    nStat_ = (uint)nStat;
    // Initialize statistic vector
    stat_.clear();
    if ( augmented_ ) {
      stat_.resize(nStat_,stat);
    }
  }

  RiskVector( const Teuchos::RCP<Vector<Real> > &vec,
              const bool augmented = false )
    : vec_(vec), augmented_(augmented), nStat_((augmented ? 1 : 0)) {
    stat_.clear();
    dual_vec1_ = vec->dual().clone();
    if (augmented) {
      stat_.resize(nStat_,0);
    }
  }
 
  RiskVector( const Teuchos::RCP<Vector<Real> > &vec,
              const std::vector<Real> &stat,
              const bool augmented = true )
    : stat_(stat), vec_(vec), augmented_(augmented), nStat_(stat.size()) {
    dual_vec1_ = vec->dual().clone();
  }

  void plus( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    if (augmented_) {
      for ( uint i = 0; i < nStat_; i++ ) {
        stat_[i] += xs.getStatistic(i);
      }
    }
    vec_->plus(*(xs.getVector()));
  }

  void scale( const Real alpha ) {
    if (augmented_) {
      for ( uint i = 0; i < nStat_; i++ ) {
        stat_[i] *= alpha;
      }
    }
    vec_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    if (augmented_) {
      for ( uint i = 0; i < nStat_; i++ ) {
        stat_[i] += alpha*xs.getStatistic(i);
      }
    }
    vec_->axpy(alpha,*(xs.getVector()));
  }

  Real dot( const Vector<Real> &x ) const {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    Real xprod = vec_->dot(*(xs.getVector()));
    if (augmented_) {
      for ( uint i = 0; i < nStat_; i++ ) {
        xprod += stat_[i]*xs.getStatistic(i);
      }
    }
    return xprod;
  }

  Real norm() const {
    return sqrt( this->dot(*this) );
  }

  const Real getStatistic(const int i = 0) const {
    TEUCHOS_TEST_FOR_EXCEPTION((i < 0 || i > (int)nStat_-1),std::invalid_argument,
      ">>> ERROR (ROL::RiskVector): index out-of-bounds in getStatistic!");
    return stat_[i];
  }

  const void getStatistic(std::vector<Real> &stat) const {
    stat.clear();
    stat.assign(stat_.begin(),stat_.end());
  }

  Teuchos::RCP<const Vector<Real> > getVector() const {
    return vec_;
  }

  Teuchos::RCP<Vector<Real> > getVector() {
    return vec_;
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    std::vector<Real> stat(nStat_,0);
    Teuchos::RCP<Vector<Real> > vec = Teuchos::rcp_dynamic_cast<Vector<Real> >(
      Teuchos::rcp_const_cast<Vector<Real> >(vec_->clone()));
    return Teuchos::rcp( new RiskVector( vec, stat, augmented_ ) );
  }

  const Vector<Real> &dual(void) const {
    dual_vec1_->set(vec_->dual());
    dual_vec_ = Teuchos::rcp(new RiskVector<Real>(dual_vec1_,stat_,augmented_));
    return *dual_vec_;
  }

  void setStatistic(const Real stat) {
    stat_.assign(nStat_,stat);
  }
 
  void setStatistic(const std::vector<Real> &stat) {
    stat_.assign(stat.begin(),stat.end());
  }
 
  void setVector(const Vector<Real>& vec) {
    vec_->set(vec);
  }
};

}

#endif
