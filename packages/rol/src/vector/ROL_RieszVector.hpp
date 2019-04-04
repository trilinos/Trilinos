
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

#ifndef ROL_RIESZVECTOR_H
#define ROL_RIESZVECTOR_H

#include <ostream>

#include "ROL_ElementwiseVector.hpp"
#include "ROL_LinearOperator.hpp"

/*
   \class ROL::RieszPrimalVector
   \brief Abstract implementation of a primal vector corresponding to
          an inner-product that involves the application of a linear 
          operator

   \class ROL::RieszDualVector
   \brief Abstract implementation of a dual vector corresponding to
          an inner-product that involves the application of a linear 
          operator inverse

*/


namespace ROL {

template<class Real>
class RieszPrimalVector;

template<class Real>
class RieszDualVector;


template <class Real>
class RieszPrimalVector : public ElementwiseVector<Real> {

  using V = Vector<Real>;
  using DualVector = RieszDualVector<Real>;
  using OP = LinearOperator<Real>;

private:

  const   ROL::Ptr<V>          v_;
  mutable ROL::Ptr<DualVector> dual_;
  const   ROL::Ptr<OP>         op_;
  mutable Real            tol_;

  mutable bool isDualInitialized_;

  void initialize_dual( void ) const {

    dual_ = ROL::makePtr<DualVector>(v_->clone(),op_,tol_);   
    op_->apply(*(dual_->getVector()),*v_,tol_);
    isDualInitialized_ = true;
  }

public:

  RieszPrimalVector( const ROL::Ptr<V>  &v, 
                     const ROL::Ptr<OP> &op,
                     Real tol=std::sqrt(ROL_EPSILON<Real>()) ) : 
    v_(v), op_(op), tol_(tol), isDualInitialized_(false) {  
  }   

  virtual ~RieszPrimalVector() {}
 
  virtual Real dot( const V &x ) const {
    if( !isDualInitialized_ ) {
      initialize_dual();
    }

    const RieszPrimalVector &ex = dynamic_cast<const RieszPrimalVector&>(x);
    return dual_->getVector()->dot(*(ex.getVector()));    
  }

  virtual ROL::Ptr<V> clone() const {
    return ROL::makePtr<RieszPrimalVector>( v_->clone(), op_, tol_ );   
  }

  virtual const V & dual() const {
    if( !isDualInitialized_ ) {
      initialize_dual();
    }
    return *dual_;
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    v_->applyUnary(f);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const V &x ) {
    const RieszPrimalVector &ex = dynamic_cast<const RieszPrimalVector&>(x);
    v_->applyBinary(f,*(ex.getVector()));
  } 

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    return v_->reduce(r); 
  }

  void setScalar( const Real C ) {
    v_->setScalar(C);
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    v_->randomize(l,u);
  }

  ROL::Ptr<V> getVector( void ) { 
    return v_;
  }

  ROL::Ptr<const V> getVector( void ) const {
    return v_;
  }
  
}; // class RieszPrimalVector



template<class Real>
class RieszDualVector : public ElementwiseVector<Real> {

  using V = Vector<Real>;
  using PrimalVector = RieszPrimalVector<Real>;
  using OP = LinearOperator<Real>;

private:

  const   ROL::Ptr<V>            v_;
  mutable ROL::Ptr<PrimalVector>  primal_;
  const   ROL::Ptr<OP>            op_;
  mutable Real               tol_;

  mutable bool isPrimalInitialized_;

  void initialize_primal( void ) const {

    primal_ = ROL::makePtr<PrimalVector>(v_->clone(),op_,tol_);   
    op_->applyInverse(*(primal_->getVector()),*v_,tol_);
    isPrimalInitialized_ = true;
  }

public:

  RieszDualVector( const ROL::Ptr<V>  &v, 
                   const ROL::Ptr<OP> &op,
                   Real tol=std::sqrt(ROL_EPSILON<Real>()) ) : 
    v_(v), op_(op), tol_(tol), isPrimalInitialized_(false) {  
  }   

  virtual ~RieszDualVector() {}
 
  virtual Real dot( const V &x ) const {
    if( !isPrimalInitialized_ ) {
      initialize_primal();
    }
   
 const RieszDualVector &ex = dynamic_cast<const RieszDualVector&>(x);
    return primal_->getVector()->dot(*(ex.getVector()));    
  }

  virtual ROL::Ptr<V> clone() const {
    return ROL::makePtr<RieszDualVector>( v_->clone(), op_, tol_ );   
  }

  virtual const V & dual() const {
    if( !isPrimalInitialized_ ) {
      initialize_primal();
    }
    return *primal_;
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    v_->applyUnary(f);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const V &x ) {
    const RieszDualVector &ex = dynamic_cast<const RieszDualVector&>(x);
    v_->applyBinary(f,*(ex.getVector()));
  } 

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    return v_->reduce(r); 
  }

  void setScalar( const Real C ) {
    v_->setScalar(C);
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    v_->randomize(l,u);
  }

  ROL::Ptr<V> getVector( void ) { 
    return v_;
  }

  ROL::Ptr<const V> getVector( void ) const {
    return v_;
  }
  
}; // class RieszDualVector



} // namespace ROL

#endif // ROL_RIESZVECTOR_H
