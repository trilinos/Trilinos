// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_COMPLEXSTDVECTOR_HPP
#define ROL_COMPLEXSTDVECTOR_HPP

#include "ROL_StdVector.hpp"
#include <algorithm>
#include <complex>
#include <functional>
#include <initializer_list>
#include <random>


namespace ROL {

template<typename Real>
inline StdVector<Real,std::complex<Real>>& complex_cast( Vector<Real>& x ) {
  return static_cast<StdVector<Real,std::complex<Real>>&>(x);
} 

template<typename Real>
inline const StdVector<Real,std::complex<Real>>& complex_cast( const Vector<Real>& x ) {
  return static_cast<const StdVector<Real,std::complex<Real>>&>(x);
} 

/** \class ROL_ComplexStdVector.hpp 
    \brief Partial specialization of ROL_StdVector for complex-valued elements */

template<typename Real>
class StdVector<Real,std::complex<Real>> : public Vector<Real> {
public:
  
  using value_type = std::complex<Real>;
  using size_type  = typename std::vector<value_type>;

  StdVector( const int dim, const value_type value=value_type{}, bool with_dual=true ) : 
    vec_(dim,value) {
  }

  StdVector( std::initializer_list<value_type> ilist ) : vec_(ilist) {
  }


  value_type& operator[] ( int i ) { return vec_[i]; }
  const value_type& operator[] ( int i ) const { return vec_[i]; }

  void set( const Vector<Real>& x ) {
    const auto& xd = complex_cast( x );
    std::copy( xd.vec_.begin(), xd.vec_.end(), vec_.begin() );
  }

  void plus( const Vector<Real>& x ) {
    const auto& xd = complex_cast( x );
    auto ptr = vec_.begin();
    auto x_ptr = xd.vec_.begin();
    while( ptr != vec_.end() ) {
      *ptr += *x_ptr;
      ++ptr; ++x_ptr;
    } 
  }

  void axpy( const Real alpha, const Vector<Real>& x ) {
    const auto& xd = complex_cast( x );
    auto ptr = vec_.begin();
    auto x_ptr = xd.vec_.begin();
    while( ptr != vec_.end() ) {
      *ptr += alpha * (*x_ptr);
      ++ptr; ++x_ptr;
    } 
  }

  void scale( const Real alpha ) { for( auto& e : vec_ ) e*=alpha; }

  Real dot( const Vector<Real>& x ) const {
    const auto& xd = complex_cast( x );
    Real result = 0;
    auto ptr = vec_.begin();
    auto x_ptr = xd.vec_.begin();
    while( ptr != vec_.end() ) {
      result += std::real( *ptr * std::conj(*x_ptr) );
      ++ptr; ++x_ptr;
    } 
    return 2*result;
  }

  Real norm() const { return std::sqrt( dot(*this) ); }

  Ptr<Vector<Real>> clone() const { return makePtr<StdVector>( dimension() ); } 

  Ptr<Vector<Real>> basis( const int i ) const { 
    auto b = makePtr<StdVector>( dimension() );
    (*b)[i] = value_type(1,0);
    return b;
  }

  int dimension() const { return static_cast<int>(vec_.size()); }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(l,u);
    for( auto& e : vec_ ) e = std::complex<Real>(dist(gen),dist(gen));
  }

  void print( std::ostream& os ) const { 
    for( auto e : vec_ ) os << e << std::endl;
    os << std::endl;
  }

private:
  std::vector<value_type> vec_;
};


template<typename Real>
using ComplexStdVector = StdVector<Real,std::complex<Real>>;


} // namespace ROL

#endif //ROL_COMPLEXSTDVECTOR_HPP

