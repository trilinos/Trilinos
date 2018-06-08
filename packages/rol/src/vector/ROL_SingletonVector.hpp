
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

  SingletonVector(const Real& value=0) : value_(value) {}

  Real getValue() const { return value_; }
  void setValue( const Real& v ) { value_=v; }

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

  void setScalar( const Real &C ) {
    value_ = C;
  }

  void print( std::ostream& os ) const {
    os << value_ << std::endl;
  }

};


} // namespace ROL




#endif // ROL_SINGLETONVECTOR_H

