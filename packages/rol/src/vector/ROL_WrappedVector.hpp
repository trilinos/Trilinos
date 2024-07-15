// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  virtual Real apply( const V& x ) const override { return vec_->apply(x); }
  
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

  virtual void setScalar( const Real C ) override { vec_->setScalar(C); }

  virtual void randomize( const Real l=0.0, const Real u=1.0 ) override {
    vec_->randomize(l,u);
  }

  virtual void print( std::ostream& os ) override { vec_->print(os); }

  const Ptr<V>& getVector() { return vec_; }
  const Ptr<const V>& getVector() const { return vec_; }

  virtual void setVector( const Ptr<const V>& vec ) const { vec_ = vec; }

}; // class PrimalScaledVector




} // namespace ROL


#endif 
