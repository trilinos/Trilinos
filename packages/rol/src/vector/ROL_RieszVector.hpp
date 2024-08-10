
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
