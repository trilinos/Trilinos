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

#ifndef ROL_CVARVECTOR_HPP
#define ROL_CVARVECTOR_HPP

#include "ROL_Vector.hpp"

namespace ROL {

template<class Real> 
class CVaRVector : public Vector<Real> {
protected:
  Real var_;
  Teuchos::RCP<Vector<Real> > vec_;

public:
  CVaRVector( Real &var, Teuchos::RCP<Vector<Real> > &vec ) : var_(var), vec_(vec) {}
  
  void plus( const Vector<Real> &x ) {
    const CVaRVector<Real> &xs = Teuchos::dyn_cast<const CVaRVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    this->var_ += xs.getVaR();
    this->vec_->plus(*(xs.getVector()));
  }   

  void scale( const Real alpha ) {
    this->var_ *= alpha;
    this->vec_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const CVaRVector<Real> &xs = Teuchos::dyn_cast<const CVaRVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    this->var_ += alpha*xs.getVaR();
    this->vec_->axpy(alpha,*(xs.getVector()));
  }

  Real dot( const Vector<Real> &x ) const {
    const CVaRVector<Real> &xs = Teuchos::dyn_cast<const CVaRVector<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    return this->var_ * xs.getVaR() + this->vec_->dot(*(xs.getVector()));
  }

  Real norm() const {
    return sqrt( this->dot(*this) );
  } 

  const Real getVaR() const { 
    return this->var_; 
  }

  Teuchos::RCP<const Vector<Real> > getVector() const { 
    return this->vec_; 
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    Real var = 0.0;
    Teuchos::RCP<Vector<Real> > vec = Teuchos::rcp_dynamic_cast<Vector<Real> >(
      Teuchos::rcp_const_cast<Vector<Real> >(this->vec_->clone()));
    return Teuchos::rcp( new CVaRVector( var, vec ) );  
  }

  void setVaR(const Real var) { 
    this->var_ = var; 
  }
  
  void setVector(const Vector<Real>& vec) { 
    this->vec_->set(vec); 
  }
};

}

#endif
