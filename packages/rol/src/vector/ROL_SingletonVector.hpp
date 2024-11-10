
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SINGLETONVECTOR_H
#define ROL_SINGLETONVECTOR_H

#include "ROL_Vector.hpp"

/** \class ROL::StdVector
    \brief Provides the ROL::Vector interface for scalar values, to be used,
           for example, with scalar constraints.
*/


namespace ROL {

template<class Real>
class SingletonVector : public Vector<Real> {

  using V  = Vector<Real>;

private:

  Real value_;

  Real getValueX( const V& x ) const { 
    return dynamic_cast<const SingletonVector<Real>&>(x).getValue(); 
  }

public:

  SingletonVector(Real value = Real(0)) : value_(value) {}

  Real getValue() const { return value_; }
  void setValue( Real v ) { value_ = v; }

  void set( const V& x ) { 
    value_ = getValueX(x);
  }

  void plus( const V& x ) {
    value_ += getValueX(x);
  }

  void axpy( const Real alpha, const V& x ) {
     value_ += alpha*getValueX(x);
  }

  void scale( const Real alpha ) {
     value_ *= alpha;
  }

  Real dot( const V& x ) const {
    Real xv = getValueX(x);
    xv *= value_;
    return xv;
  }

  Real norm() const {
    return std::abs(value_);
  }
  
  ROL::Ptr<V> clone() const {
    return ROL::makePtr<SingletonVector>(0);
  }
  
  ROL::Ptr<V> basis(const int i) const {
    ROL_TEST_FOR_EXCEPTION( i >= 1 || i < 0,
                                std::invalid_argument,
                                "Error: Basis index must be between 0 and vector dimension." );
    return ROL::makePtr<SingletonVector>(1);
  }

  int dimension() const { return 1; };

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    value_ = f.apply(value_);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const V& x ) {
    value_ = f.apply(value_,getValueX(x));
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    return value_;
  }

  void setScalar( const Real C ) {
    value_ = C;
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    Real a = (u-l);
    Real b = l;
    Real x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    value_ = a*x + b;
  }

  void print( std::ostream& os ) const {
    os << value_ << std::endl;
  }

};


} // namespace ROL




#endif // ROL_SINGLETONVECTOR_H

