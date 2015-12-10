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
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real> 
class RiskVector : public Vector<Real> {
protected:
  Real stat_;
  Teuchos::RCP<Vector<Real> > vec_;
  bool augmented_;

  mutable Teuchos::RCP<Vector<Real> > dual_vec1_;
  mutable Teuchos::RCP<RiskVector<Real> > dual_vec_;

public:
  RiskVector( Teuchos::ParameterList &parlist,
        const Teuchos::RCP<Vector<Real> > &vec,
        const Real stat = 1. )
    : stat_(stat), vec_(vec), augmented_(false) {
    dual_vec1_ = vec->dual().clone();
    std::string type = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
    if ( type == "CVaR" || type == "HMCR" ||
         type == "Moreau-Yosida CVaR" ||
         type == "Log-Exponential Quadrangle" ||
         type == "Log-Quantile Quadrangle" ||
         type == "Quantile-Based Quadrangle" ||
         type == "Truncated Mean Quadrangle" ) {
      augmented_ = true;
    }
  }

  RiskVector( const Teuchos::RCP<Vector<Real> > &vec,
              const bool augmented = false,
              const Real stat = 0.)
    : stat_(stat), vec_(vec), augmented_(augmented) {
    dual_vec1_ = vec->dual().clone();
  }
  
  void plus( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    if (augmented_) { stat_ += xs.getStatistic(); }
    vec_->plus(*(xs.getVector()));
  }   

  void scale( const Real alpha ) {
    if (augmented_) { stat_ *= alpha; }
    this->vec_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    if (augmented_) { stat_ += alpha*xs.getStatistic(); }
    vec_->axpy(alpha,*(xs.getVector()));
  }

  Real dot( const Vector<Real> &x ) const {
    const RiskVector<Real> &xs = Teuchos::dyn_cast<const RiskVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    Real xprod = vec_->dot(*(xs.getVector()));
    if (augmented_) { xprod += stat_*xs.getStatistic(); }
    return xprod;
  }

  Real norm() const {
    return sqrt( this->dot(*this) );
  } 

  const Real getStatistic() const { 
    return stat_; 
  }

  Teuchos::RCP<const Vector<Real> > getVector() const { 
    return vec_; 
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    Real stat = 0.0;
    Teuchos::RCP<Vector<Real> > vec = Teuchos::rcp_dynamic_cast<Vector<Real> >(
      Teuchos::rcp_const_cast<Vector<Real> >(vec_->clone()));
    return Teuchos::rcp( new RiskVector( vec, augmented_, stat ) );  
  }

  const Vector<Real> &dual(void) const {
    dual_vec1_->set(vec_->dual());
    dual_vec_ = Teuchos::rcp(new RiskVector<Real>(dual_vec1_,augmented_,stat_));
    return *dual_vec_;
  }

  void setStatistic(const Real stat) { 
    stat_ = stat; 
  }
  
  void setVector(const Vector<Real>& vec) { 
    vec_->set(vec); 
  }
};

}

#endif
