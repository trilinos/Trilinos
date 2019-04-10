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

#ifndef ROL_STDARRAY_H
#define ROL_STDARRAY_H

#include <algorithm>
#include <array>
#include <utility>
#include <random>
#include "ROL_Vector.hpp"

/** \class ROL::StdArray
    \brief Provides the std::array implementation of the ROL::Vector interface.
*/


namespace ROL {

template<typename Real, std::size_t N>
class StdArray : public Vector<Real> {
public:
  
  StdArray( std::array<Real,N>&& array_in  ) : data(std::move(array_in)) {}
  StdArray() = default;


  inline Real& operator[] ( std::size_t i ) { return data[i]; }
  inline const Real& operator[] ( std::size_t i ) const { return data[i]; }

  std::array<Real,N>& getArray() { return data; }
  const std::array<Real,N>& getArray() const { return data; }

  void set( const Vector<Real> &x ) {
    auto& ex = static_cast<const StdArray&>(x).getArray();    
    std::copy(ex.begin(),ex.end(),data.begin());
  }

  void plus( const Vector<Real> &x ) {
    auto& ex = static_cast<const StdArray&>(x).getArray();    
    std::transform(ex.begin(),ex.end(),data.begin(),data.begin(),std::plus<Real>{});
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    auto& ex = static_cast<const StdArray&>(x).getArray();    
    std::transform(ex.begin(),ex.end(),data.begin(),data.begin(),[alpha](Real x, Real y){ return alpha*x+y; });
  }

  void scale( const Real alpha ) {
    for( auto& e : data ) e *= alpha;
  }

  virtual Real dot( const Vector<Real> &x ) const {
    Real result = 0;
    auto& ex = static_cast<const StdArray&>(x).getArray();    
    std::inner_product(ex.begin(),ex.end(),data.begin(),result);
    return result;
  }

  Real norm() const {
    Real norm_squared = 0;
    for( auto e: data ) norm_squared += (e*e);
    return std::sqrt(norm_squared);
  }

  virtual Ptr<Vector<Real>> clone() const {
    return makePtr<StdArray>(std::move(std::array<Real,N>()));  
  }

  Ptr<Vector<Real>> basis( const int i ) const {
    std::array<Real,N> b_array;
    b_array.fill(0);
    b_array[i] = Real(1);
    return makePtr<StdArray>(std::move(b_array));
  }

  int dimension() const { return N; }

  void zero() { data.fill(0); }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    for( auto& e : data ) e = f.apply(e);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, 
                    const Vector<Real> &x ) {
    auto& ex = static_cast<const StdArray&>(x).getArray();    
    std::transform(ex.begin(),ex.end(),data.begin(),data.begin(),
                   [&f](Real a, Real b){ return f.apply(a,b);});
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    for( auto e: data ) r.reduce(e,result);
    return result;
  }

  void setScalar( const Real alpha ) { data.fill(alpha); }

  void randomize( const Real l = -1.0, const Real u = 1.0 ) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<Real> dis(l, u);
    for( auto& e : data ) e = dis(gen);
  }

  virtual void print( std::ostream &outStream ) const {
    for( auto e: data ) outStream << e << " ";
    outStream << std::endl;
  }

private:

  const std::array<Real,N> _array( const Vector<Real>& x ) {
    return static_cast<const StdArray&>(x).getArray();    
  }

  std::array<Real,N> data;

}; // class StdArray


} // namespace ROL

#endif
