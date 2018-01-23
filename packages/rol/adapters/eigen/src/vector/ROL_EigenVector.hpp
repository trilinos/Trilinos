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


#pragma once

#include "ROL_ElementwiseVector.hpp"
#include <Eigen/Core>

/** \class ROL::EigenVector
    \brief Provides the Eigen 3 vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template<class Real>
class Eigen3Vector : public ElementwiseVector<Real> {

  using V = Vector<Real>;

  using EV = Eigen::Matrix<Real,Eigen::Dynamic,1>;

  using UF = Elementwise::UnaryFunction<Real>;
  using BF = Elementwise::BinaryFunction<Real>;
  using RO = Elementwise::ReductionOp<Real>;

private:

  ROL::Ptr<EV> vec_;

  int dim_;

public:

  Eigen3Vector( const ROL::Ptr<EV> &vec ) : vec_(vec), dim_(vec->size()) {
  }
 
  Eigen3Vector( int dim, bool zeroOut=false ) : dim_(dim) {
    if( zeroOut ) {
      vec_ = ROL::makePtr<EV>(dim_);
    }
    else {
      vec_ = ROL::makePtr<EV>(dim_);
    }
  }

  void applyUnary( const UF &f ) {
    for( int i=0; i<dim_; ++i ) 
      (*vec_)(i) = f.apply((*vec_)(i));
  }

  void applyBinary( const BF &f, const V &x ) {
    auto ex = dynamic_cast<const Eigen3Vector&>(x);
    for( int i=0; i<dim_; ++i ) 
      (*vec_)(i) = f.apply((*vec_)(i),ex(i));
  }

  Real reduce( const RO &r ) const {
    Real result = r.initialValue();
    for( int i=0; i<dim_; ++i ) 
      r.reduce((*vec_)(i),result);
    return result;
  }

  int dimension() const {
    return dim_;
  }

  ROL::Ptr<V> basis( const int i ) const {
    auto data = ROL::makePtr<EV>(dim_,true);
    (*data)(i) = static_cast<Real>(1.0);
    return ROL::makePtr<Eigen3Vector>(data);
  }

  ROL::Ptr<V> clone() const {
    return ROL::makePtr<Eigen3Vector>(dim_); 
  }

  ROL::Ptr<EV> getVector() {
    return vec_;
  }

  ROL::Ptr<const EV> getVector() const {
    return vec_;
  }

  Real& operator() ( int i ) {
    return (*vec_)(i);
  }

  const Real& operator() ( int i ) const {
    return (*vec_)(i);
  }

  void print( std::ostream &os ) const {
    os << *vec_ << std::endl;
  }
  
}; // class Vector

} // namespace ROL

