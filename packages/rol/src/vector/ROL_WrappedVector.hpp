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
#ifndef ROL_WRAPPEDVECTOR_HPP
#define ROL_WRAPPEDVECTOR_HPP


namespace ROL {

/** @ingroup la_group 
    \class ROL::WrappedVector
    \brief Provides an interface layer which encapulates a pointer to a 
           ROL::Vector and has the default behavior of calling its
           member Ptr<Vector> object. Makes creating derived classes
           with this idiom simpler as they need only override the 
           methods where the desired implementation differs from
           the member Ptr<Vector>. For example, vectors which 
           have a diagonal scaling that defines their inner product and
           dual spaces can derive from this class need overload only
           the methods basis, clone, dual, and dot.

*/

template<typename Real>
class PrimalScaledVector {

  using V     = Vector<Real>;
  using VPrim = PrimalScaledVector<Real>;
  using VDual = DualScaledVector<Real>;

private:

  mutable Ptr<V>   vec_;


public: 

  WrappedVector( const Ptr<V>& vec ) : vec_(vec) {}

  virtual ~WrappedVector() {}

  virtual void plus( const V& x ) override { vec_->plus(x); }
  virtual void scale( const Real alpha ) override { vec_->scale(alpha); }
  
  virtual Real dot( const V& x ) const override {  return vec_->dot(x); }

  virtual Real norm() const override { return std::sqrt( this->dot(*this) ); }

  virtual Ptr<V> clone() const override { 
    return makePtr<WrappedVector>( vec_->clone() );
  }

  virtual void axpy( const Real alpha, const V& x ) override {
    vec_->axpy( alpha, x );
  }

  virtual Ptr<V> basis( const int i ) const override {
    return makePtr<VPrim>( vec_->basis(i) );
  }

  virtual int dimension() const override { return vec_->dimension(); }

  virtual void set( const V& x ) override { vec_->set(x); }

  virtual void const V& dual() const override { return vec_->dual(); }
  
  virtual void applyUnary( const Elementwise::UnaryFunction<Real>& f ) override {
    vec_->applyUnary(f);
  }

  virtual void applyBinary( const Elementwise::BinaryFunction<Real>& f, 
                            const V& x ) override {
    vec_->applyBinary(f,x);
  }
  
  virtual Real reduce( const Elementwise::ReductionOp<Real>& r ) const override {
    return vec_->reduce(r);
  }

  virtual void setScalar( const Real& c ) override { vec_->setScalar(c); }

  virtual void print( std::ostream& os ) override { vec_->print(os); }

  const Ptr<V>& getVector() { return vec_; }
  const Ptr<const V>& getVector() const { return vec_; }

  virtual void setVector( const Ptr<const V>& vec ) const { vec_ = vec; }

}; // class PrimalScaledVector




} // namespace ROL


#endif 
