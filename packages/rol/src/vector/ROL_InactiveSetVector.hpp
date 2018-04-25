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

#ifndef ROL_INACTIVE_SET_VECTOR_HPP
#define ROL_INACTIVE_SET_VECTOR_HPP

#include "ROL_Vector.hpp"
#include "ROL_VectorWorkspace.hpp"
#include "ROL_BoundConstraint.hpp"


/** @ingroup la_group
    \class ROL::InActiveVector
    \brief Defines the a Vector which has a diagonally scaled dot 
           product that neglects active set elements

           Used to simplify Semi-smooth Newton method implementation
*/

namespace ROL {

template<typename Real> class InactiveSet_PrimalVector;
template<typename Real> class InactiveSet_DualVector;


template<typename Real>
class InactiveSet_PrimalVector : public Vector<Real> {

  using V      = Vector<Real>;
  using Primal = InactiveSet_PrimalVector<Real>;
  using Dual   = InactiveSet_DualVector<Real>;
  using Bnd    = BoundConstraint<Real>;

private:

  Ptr<V>          vec_;
  Ptr<V>          scaling_vec_;
  Ptr<Bnd>        bnd_;
  
  VectorWorkspace<Real> workspace_;  

  Elementwise::Multiply<Real> mult_;

public: 

  InactiveSet_PrimalVector( const Ptr<V>& vec,
                            const Ptr<V>& scaling_vec, 
                            const Ptr<Bnd>& bnd ) :
    vec_(vec), scaling_vec_(scaling_vec), bnd_(bnd) {}

  virtual ~InactiveSet_PrimalVector() {}

  void plus( const V& x ) override { vec_->plus(x); }
  void scale( const Real alpha ) override { vec_->scale(alpha); }
  


  Real dot( const V& x ) const override {

    auto y = workspace_.copy(x);

    y->applyBinary( mult_, *scaling_ );
    bnd_->pruneActive( *y, *vec_ );

    return y->dot(*vec_);    
  } 

  Real norm() const override { return std::sqrt(this->dot(*this)); }

  Ptr<V> clone() const override {
    return makePtr<Primal>( vec_->clone(), scaling_vec_, bnd_ );
  }

  void axpy( const Real alpha, const Vector& x ) override { 
    vec_->axpy( alpha, x );
  }

  Ptr<V> basis( const int i ) const override { 
    return makePtr<Primal>( vec_->basis(i), scaling_vec_, bnd_ ); 
  }

  int dimension() const override { return vec_->dimension(); }
 
  void set( const V& x ) override { vec_->set(x); }

  void const V& dual() const override {
    auto dual_vec = workspace_.copy( *vec_ );  
    dual_vec->applyBinary( mult_, *scaling_vec_ );
    return makePtr<Dual>( dual_vec, scaling_vec_, bnd_ );
  } 

  void applyUnary( const Elementwise::UnaryFunction<Real>& f ) override {
    vec_->applyUnary( f );  
  }

  void applyBinary( const Elementwise::UnaryFunction<Real>& f, const V& x ) override {
    vec_->applyBinary( f, x );
  }

  Real reduce( const Elementwise::ReductionOp<Real>& r ) const override {
    return vec_->reduce( r );
  }

  void setScalar( const Real& c ) override { vec_->set(c); }

  void print( std::ostream& os ) const override { vec_->print(os); }

}; // class InactiveSet_PrimalVector 





template<typename Real>
class InactiveSet_DualVector : public Vector<Real> {

  using V      = Vector<Real>;
  using Primal = InactiveSet_PrimalVector<Real>;
  using Dual   = InactiveSet_DualVector<Real>;
  using Bnd    = BoundConstraint<Real>;

private:

  Ptr<V>          vec_;
  Ptr<V>          scaling_vec_;
  Ptr<Bnd>        bnd_;
  
  VectorWorkspace<Real> workspace_;  

  Elementwise::Divide<Real>   div_;

public: 

  InactiveSet_DualVector( const Ptr<V>& vec,
                          const Ptr<V>& scaling_vec, 
                          const Ptr<Bnd>& bnd ) :
    vec_(vec), scaling_vec_(scaling_vec), bnd_(bnd) {}

  virtual ~InactiveSet_PrimalVector() {}

  void plus( const V& x ) override { vec_->plus(x); }
  void scale( const Real alpha ) override { vec_->scale(alpha); }
  
  Real dot( const V& x ) const override {

    auto y = workspace_.copy(x);

    y->applyBinary( div_, *scaling_ );
    bnd_->pruneActive( *y, *vec_ );

    return y->dot(*vec_);    
  } 

  Real norm() const override { return std::sqrt(this->dot(*this)); }

  Ptr<V> clone() const override {
    return makePtr<Primal>( vec_->clone(), scaling_vec_, bnd_ );
  }

  void axpy( const Real alpha, const Vector& x ) override { 
    vec_->axpy( alpha, x );
  }

  Ptr<V> basis( const int i ) const override { 
    return makePtr<Primal>( vec_->basis(i), scaling_vec_, bnd_ ); 
  }

  int dimension() const override { return vec_->dimension(); }
 
  void set( const V& x ) override { vec_->set(x); }

  void const V& dual() const override {
    auto primal_vec = workspace_.copy( *vec_ );  
    primal_vec->applyBinary( div_, *scaling_vec_ );
    return makePtr<Primal>( primal_vec, scaling_vec_, bnd_ );
  } 

  void applyUnary( const Elementwise::UnaryFunction<Real>& f ) override {
    vec_->applyUnary( f );  
  }

  void applyBinary( const Elementwise::UnaryFunction<Real>& f, const V& x ) override {
    vec_->applyBinary( f, x );
  }

  Real reduce( const Elementwise::ReductionOp<Real>& r ) const override {
    return vec_->reduce( r );
  }

  void setScalar( const Real& c ) override { vec_->set(c); }

  void print( std::ostream& os ) const override { vec_->print(os); }

}; // class InactiveSet_DualVector 


} // namespace ROL



#endif // ROL_INACTIVE_SET_VECTOR_HPP
