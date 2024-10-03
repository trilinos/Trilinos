// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDARRAY_H
#define ROL_STDARRAY_H

#include <algorithm>
#include <array>
#include <utility>
#include <random>
#include <numeric>
#include "ROL_Vector.hpp"

/** \class ROL::StdArray
    \brief Provides the std::array implementation of the ROL::Vector interface.
*/


namespace ROL {

template<typename Real, std::size_t array_size, std::size_t pool_size=100u>
class StdArray : public Vector<Real> {
public:

  using data_type = std::array<Real,array_size>;
  
  StdArray() { 
    for( auto& vptr : pool_ptr ) {
      if( getCount(vptr) < 2 ) {
        data = vptr;
        break;
      }
    }
    if( is_nullPtr(data) ) {
      data = makePtr<std::array<Real,array_size>>();
    } 
  }
 

  inline Real& operator[] ( std::size_t i ) { return (*data)[i]; }
  inline const Real& operator[] ( std::size_t i ) const { return (*data)[i]; }

  std::array<Real,array_size>& get_array() { return *data; }
  const std::array<Real,array_size>& get_array() const { return *data; }

  void set( const Vector<Real> &x ) {
    const auto& ex = _array(x);
    std::copy(ex.begin(),ex.end(),data->begin());
  }

  void plus( const Vector<Real> &x ) {
    const auto& ex = _array(x);
    std::transform(ex.begin(),ex.end(),data->begin(),data->begin(),std::plus<Real>{});
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const auto& ex = _array(x);
    std::transform(ex.begin(),ex.end(),data->begin(),data->begin(),[alpha](Real x, Real y){ return alpha*x+y; });
  }

  void scale( const Real alpha ) {
    for( auto& e : *data ) e *= alpha;
  }

  virtual Real dot( const Vector<Real> &x ) const {
    Real result = 0;
    const auto& ex = _array(x);
    std::inner_product(ex.begin(),ex.end(),data->begin(),result);
    return result;
  }

  Real norm() const {
    Real norm_squared = 0;
    for( auto e: *data ) norm_squared += (e*e);
    return std::sqrt(norm_squared);
  }

  virtual Ptr<Vector<Real>> clone() const {
    return makePtr<StdArray>();  
  }

  Ptr<Vector<Real>> basis( const int i ) const {
    auto  b_ptr = clone();
    auto& b_ref = static_cast<StdArray&>(*b_ptr);
    b_ref.zero();
    b_ref[i] = Real(1);
    return b_ptr;
  }

  int dimension() const { return static_cast<int>(array_size); }

  void zero() { data->fill(0); }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    for( auto& e : *data ) e = f.apply(e);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, 
                    const Vector<Real> &x ) {
    const auto& ex = _array(x);
    std::transform(ex.begin(),ex.end(),data->begin(),data->begin(),
                   [&f](Real a, Real b){ return f.apply(a,b);});
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    for( auto e: *data ) r.reduce(e,result);
    return result;
  }

  void setScalar( const Real alpha ) { data->fill(alpha); }

  void randomize( const Real l = -1.0, const Real u = 1.0 ) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<Real> dis(l, u);
    for( auto& e : *data ) e = dis(gen);
  }

  virtual void print( std::ostream &outStream ) const {
    for( auto e: *data ) outStream << e << " ";
    outStream << std::endl;
  }

  static void initialize_pool() {
    for( std::size_t i=0; i<array_size; ++i ) pool_ptr[i] = makePtrFromRef(pool[i]);
  }

  // Count how many objects in the pool are currently being used 
  static std::size_t pool_count() {
    std::size_t count = 0u;
    for( auto& vptr : pool_ptr ) count += ( getCount(vptr)>1 );
    return count;
  }

private:

  StdArray( Ptr<std::array<Real,array_size>> p ) : data(p) {}

  const std::array<Real,array_size>& _array( const Vector<Real>& x ) const {
    return static_cast<const StdArray&>(x).get_array();    
  }

  Ptr<std::array<Real,array_size>> data;

  // Allocate scratch space at compile time
  static std::array<std::array<Real,array_size>,pool_size>      pool; 
  static std::array<Ptr<std::array<Real,array_size>>,pool_size> pool_ptr;

}; // class StdArray

template<typename Real, std::size_t array_size, std::size_t pool_size>
std::array<std::array<Real,array_size>,pool_size> StdArray<Real,array_size,pool_size>::pool; 

template<typename Real, std::size_t array_size, std::size_t pool_size>
std::array<Ptr<std::array<Real,array_size>>,pool_size>  StdArray<Real,array_size,pool_size>::pool_ptr; 

} // namespace ROL

#endif
