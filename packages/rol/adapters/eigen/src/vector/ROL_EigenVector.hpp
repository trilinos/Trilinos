// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    vec_ = ROL::makePtr<EV>(dim_);
    if( zeroOut ) vec_->setZero();
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
    auto data = ROL::makePtr<EV>(dim_);
    data->setZero();
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

