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

#ifndef ROL_VECTOR_SIMOPT_HPP
#define ROL_VECTOR_SIMOPT_HPP

#include "ROL_Vector.hpp"

/** @ingroup la_group
    \class ROL::Vector_SimOpt
    \brief Defines the linear algebra or vector space interface for simulation-based optimization.
*/

namespace ROL {

template<class Real> 
class Vector_SimOpt : public Vector<Real> {
private:
  Teuchos::RCP<Vector<Real> > vec1_;
  Teuchos::RCP<Vector<Real> > vec2_;

public:
  Vector_SimOpt( const Teuchos::RCP<Vector<Real> > &vec1, const Teuchos::RCP<Vector<Real> > &vec2 ) 
    : vec1_(vec1), vec2_(vec2) {}
  
  void plus( const Vector<Real> &x ) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    vec1_->plus(*(xs.get_1()));
    vec2_->plus(*(xs.get_2()));
  }   

  void scale( const Real alpha ) {
    vec1_->scale(alpha);
    vec2_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    vec1_->axpy(alpha,*(xs.get_1()));
    vec2_->axpy(alpha,*(xs.get_2()));
  }

  Real dot( const Vector<Real> &x ) const {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    return vec1_->dot(*(xs.get_1())) + vec2_->dot(*(xs.get_2()));
  }

  Real norm() const {
    Real norm1 = vec1_->norm();
    Real norm2 = vec2_->norm();
    return sqrt( norm1*norm1 + norm2*norm2 );
  } 

  Teuchos::RCP<Vector<Real> > clone() const {
    return Teuchos::rcp( new Vector_SimOpt(vec1_->clone(),vec2_->clone()) );  
  }

  Teuchos::RCP<Vector<Real> > basis( const int i )  const {
    int n1 = (vec1_)->dimension();
    if ( i < n1 ) {
      Teuchos::RCP<Vector<Real> > e1 = (vec1_)->basis(i);
      Teuchos::RCP<Vector<Real> > e2 = (vec2_)->clone(); e2->zero();
      Teuchos::RCP<Vector<Real> > e  = Teuchos::rcp(new Vector_SimOpt<Real>(e1,e2));
      return e;
    }
    else {
      Teuchos::RCP<Vector<Real> > e1 = (vec1_)->clone(); e1->zero();
      Teuchos::RCP<Vector<Real> > e2 = (vec2_)->basis(i-n1);
      Teuchos::RCP<Vector<Real> > e  = Teuchos::rcp(new Vector_SimOpt<Real>(e1,e2));
      return e;
    }
  }

  int dimension() const {
    return (vec1_)->dimension() + (vec2_)->dimension();
  }

  Teuchos::RCP<const Vector<Real> > get_1() const { 
    return vec1_; 
  }

  Teuchos::RCP<const Vector<Real> > get_2() const { 
    return vec2_; 
  }

  void set_1(const Vector<Real>& vec) { 
    vec1_->set(vec);
  }
  
  void set_2(const Vector<Real>& vec) { 
    vec2_->set(vec); 
  }
};

}

#endif
