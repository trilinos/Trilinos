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

#include "ROL_StdVector.hpp"
#include "ROL_RiskMeasureInfo.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real> 
class RiskVector : public Vector<Real> {
private:
  Teuchos::RCP<std::vector<Real> > stat_;
  Teuchos::RCP<StdVector<Real> >   stat_vec_;
  Teuchos::RCP<Vector<Real> >      vec_;

  bool augmented_;
  int nStat_;

  mutable bool isDualInitialized_;
  mutable Teuchos::RCP<Vector<Real> > dual_vec1_;
  mutable Teuchos::RCP<RiskVector<Real> > dual_vec_;

public:
  RiskVector( Teuchos::ParameterList &parlist,
        const Teuchos::RCP<Vector<Real> > &vec,
        const Real stat = 1 )
    : stat_(Teuchos::null), stat_vec_(Teuchos::null), vec_(vec),
      augmented_(false), nStat_(0), isDualInitialized_(false) {
    // Get risk measure information
    std::string name;
    std::vector<Real> lower, upper;
    bool activated(false);
    int nStat(0);
    RiskMeasureInfo<Real>(parlist,name,nStat,lower,upper,activated);
    augmented_ = (nStat > 0) ? true : false;
    nStat_ = nStat;
    // Initialize statistic vector
    if (augmented_) {
      stat_ = Teuchos::rcp(new std::vector<Real>(nStat_,stat));
      stat_vec_ = Teuchos::rcp(new StdVector<Real>(stat_));
    }
  }

  RiskVector( const Teuchos::RCP<Vector<Real> > &vec,
              const bool augmented = false )
    : stat_(Teuchos::null), stat_vec_(Teuchos::null), vec_(vec),
      augmented_(augmented), nStat_((augmented ? 1 : 0)), isDualInitialized_(false) {
    if (augmented) {
      stat_ = Teuchos::rcp(new std::vector<Real>(nStat_,0));
      stat_vec_ = Teuchos::rcp(new StdVector<Real>(stat_));
    }
  }
 
  RiskVector( const Teuchos::RCP<Vector<Real> > &vec,
              const std::vector<Real> &stat,
              const bool augmented = true )
    : stat_(Teuchos::null), stat_vec_(Teuchos::null), vec_(vec),
      augmented_(augmented), nStat_(stat.size()), isDualInitialized_(false) {
    if (augmented) {
      stat_ = Teuchos::rcp(new std::vector<Real>(stat));
      stat_vec_ = Teuchos::rcp(new StdVector<Real>(stat_));
    }
  }

  void set( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    vec_->set(*(xs.getVector()));
    if (augmented_) {
      stat_vec_->set(*(xs.getStatistic()));
    }
  }

  void plus( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    vec_->plus(*(xs.getVector()));
    if (augmented_) {
      stat_vec_->plus(*(xs.getStatistic()));
    }
  }

  void scale( const Real alpha ) {
    vec_->scale(alpha);
    if (augmented_) {
      stat_vec_->scale(alpha);
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    vec_->axpy(alpha,*(xs.getVector()));
    if (augmented_) {
      stat_vec_->axpy(alpha,*(xs.getStatistic()));
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    Real val = vec_->dot(*(xs.getVector()));
    if (augmented_) {
      val += stat_vec_->dot(*(xs.getStatistic()));
    }
    return val;
  }

  Real norm(void) const {
    return sqrt( dot(*this) );
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    std::vector<Real> stat(nStat_,static_cast<Real>(0));
    return Teuchos::rcp(new RiskVector(vec_->clone(),stat,augmented_));
  }

  const Vector<Real> &dual(void) const {
    // Initialize dual vectors if not already initialized
    if ( !isDualInitialized_ ) {
      dual_vec1_ = vec_->dual().clone();
      std::vector<Real> stat(nStat_,static_cast<Real>(0));
      dual_vec_  = Teuchos::rcp(new RiskVector<Real>(dual_vec1_,stat,augmented_));
      isDualInitialized_ = true;
    }
    // Set vector component 
    dual_vec1_->set(vec_->dual());
    // Set statistic component
    Teuchos::rcp_dynamic_cast<RiskVector<Real> >(dual_vec_)->setStatistic(*stat_);
    // Return dual vector
    return *dual_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i )  const {
    Teuchos::RCP<Vector<Real> > e1;
    std::vector<Real> e2(nStat_,static_cast<Real>(0));
    int n1 = vec_->dimension(), n2 = stat_vec_->dimension();
    if ( i < n1 ) {
      e1 = vec_->basis(i);
    }
    else if (i >= n1 && i < n1+n2) {
      e1 = vec_->clone(); e1->zero();
      e2[i-n1] = static_cast<Real>(1);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> ERROR (ROL::RiskVector::Basis): index is out of bounds.");
    }
    return Teuchos::rcp(new RiskVector<Real>(e1,e2,augmented_));
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    vec_->applyUnary(f);
    if (augmented_) {
      stat_vec_->applyUnary(f);
    }
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    vec_->applyBinary(f,*xs.getVector());
    if (augmented_) {
      stat_vec_->applyBinary(f,*xs.getStatistic());
    }
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    r.reduce(vec_->reduce(r),result);
    if (augmented_) {
      r.reduce(stat_vec_->reduce(r),result);
    }
    return result;
  }

  int dimension(void) const {
    int dim = vec_->dimension();
    if (augmented_) {
      dim += stat_vec_->dimension();
    }
    return dim;
  }

  /***************************************************************************/
  /************ ROL VECTOR ACCESSOR FUNCTIONS ********************************/
  /***************************************************************************/
  Teuchos::RCP<const StdVector<Real> > getStatistic(void) const {
    return stat_vec_;
  }

  Teuchos::RCP<StdVector<Real> > getStatistic(void) {
    return stat_vec_;
  }

  Teuchos::RCP<const Vector<Real> > getVector(void) const {
    return vec_;
  }

  Teuchos::RCP<Vector<Real> > getVector(void) {
    return vec_;
  }

  /***************************************************************************/
  /************ COMPONENT ACCESSOR FUNCTIONS *********************************/
  /***************************************************************************/
  const Real getStatistic(const int i) const {
    TEUCHOS_TEST_FOR_EXCEPTION((i < 0 || i > (int)nStat_-1),std::invalid_argument,
      ">>> ERROR (ROL::RiskVector::getStatistic): index out-of-bounds.");
    TEUCHOS_TEST_FOR_EXCEPTION(!augmented_,std::invalid_argument,
      ">>> ERROR (ROL::RiskVector::getStatistic): vector is not augmented.");
    return (*stat_)[i];
  }

  void getStatistic(std::vector<Real> &stat) const {
    stat.clear();
    if ( augmented_ ) {
      stat.assign(stat_->begin(),stat_->end());
    }
  }

  void setStatistic(const Real stat) {
    if ( augmented_ ) {
      stat_->assign(nStat_,stat);
    }
  }
 
  void setStatistic(const std::vector<Real> &stat) {
    if ( augmented_ ) {
      stat_->assign(stat.begin(),stat.end());
    }
  }
 
  void setVector(const Vector<Real>& vec) {
    vec_->set(vec);
  }
};

}

#endif
