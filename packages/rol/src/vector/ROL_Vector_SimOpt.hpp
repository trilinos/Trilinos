// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VECTOR_SIMOPT_HPP
#define ROL_VECTOR_SIMOPT_HPP

#include "ROL_Vector.hpp"

/** @ingroup la_group
    \class ROL::Vector_SimOpt
    \brief Defines the linear algebra or vector space interface for simulation-based optimization.
*/

namespace ROL {

template<class Real> 
class Vector_SimOpt : public Vector<Real> {
private:
  ROL::Ptr<Vector<Real> > vec1_;
  ROL::Ptr<Vector<Real> > vec2_;
  mutable ROL::Ptr<Vector<Real> > dual_vec1_;
  mutable ROL::Ptr<Vector<Real> > dual_vec2_;
  mutable ROL::Ptr<Vector_SimOpt<Real> > dual_vec_;

public:
  Vector_SimOpt( const ROL::Ptr<Vector<Real> > &vec1, const ROL::Ptr<Vector<Real> > &vec2 ) 
    : vec1_(vec1), vec2_(vec2) {
    dual_vec1_ = (vec1_->dual()).clone();
    dual_vec2_ = (vec2_->dual()).clone();
  }

  void plus( const Vector<Real> &x ) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    vec1_->plus(*(xs.get_1()));
    vec2_->plus(*(xs.get_2()));
  }   

  void scale( const Real alpha ) {
    vec1_->scale(alpha);
    vec2_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    vec1_->axpy(alpha,*(xs.get_1()));
    vec2_->axpy(alpha,*(xs.get_2()));
  }

  Real dot( const Vector<Real> &x ) const {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    return vec1_->dot(*(xs.get_1())) + vec2_->dot(*(xs.get_2()));
  }

  Real norm() const {
    Real norm1 = vec1_->norm();
    Real norm2 = vec2_->norm();
    return sqrt( norm1*norm1 + norm2*norm2 );
  } 

  ROL::Ptr<Vector<Real> > clone() const {
    return ROL::makePtr<Vector_SimOpt>(vec1_->clone(),vec2_->clone());  
  }

  const Vector<Real> & dual(void) const {
    dual_vec1_->set(vec1_->dual());
    dual_vec2_->set(vec2_->dual());
    dual_vec_ = ROL::makePtr<Vector_SimOpt<Real>>(dual_vec1_,dual_vec2_); 
    return *dual_vec_;
  }

  Real apply( const Vector<Real> &x ) const {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    return vec1_->apply(*(xs.get_1())) + vec2_->apply(*(xs.get_2()));
  }

  ROL::Ptr<Vector<Real> > basis( const int i )  const {
    int n1 = (vec1_)->dimension();
    if ( i < n1 ) {
      ROL::Ptr<Vector<Real> > e1 = (vec1_)->basis(i);
      ROL::Ptr<Vector<Real> > e2 = (vec2_)->clone(); e2->zero();
      ROL::Ptr<Vector<Real> > e  = ROL::makePtr<Vector_SimOpt<Real>>(e1,e2);
      return e;
    }
    else {
      ROL::Ptr<Vector<Real> > e1 = (vec1_)->clone(); e1->zero();
      ROL::Ptr<Vector<Real> > e2 = (vec2_)->basis(i-n1);
      ROL::Ptr<Vector<Real> > e  = ROL::makePtr<Vector_SimOpt<Real>>(e1,e2);
      return e;
    }
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {

    vec1_->applyUnary(f);
    vec2_->applyUnary(f);

  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);

    vec1_->applyBinary(f,*xs.get_1());
    vec2_->applyBinary(f,*xs.get_2());
  
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {

    Real result = r.initialValue();
    r.reduce(vec1_->reduce(r),result);
    r.reduce(vec2_->reduce(r),result);
    return result;
  }

  void setScalar( const Real C ) {
    vec1_->setScalar(C);
    vec2_->setScalar(C);
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    vec1_->randomize(l,u);
    vec2_->randomize(l,u);
  }


  int dimension() const {
    return (vec1_)->dimension() + (vec2_)->dimension();
  }

  ROL::Ptr<const Vector<Real> > get_1() const { 
    return vec1_; 
  }

  ROL::Ptr<const Vector<Real> > get_2() const { 
    return vec2_; 
  }

  ROL::Ptr<Vector<Real> > get_1() { 
    return vec1_; 
  }

  ROL::Ptr<Vector<Real> > get_2() { 
    return vec2_; 
  }

  void set_1(const Vector<Real>& vec) { 
    vec1_->set(vec);
  }
  
  void set_2(const Vector<Real>& vec) { 
    vec2_->set(vec); 
  }

  void print( std::ostream &outStream ) const {
    outStream << "Sim: ";
    vec1_->print(outStream);
    outStream << "Opt: ";
    vec2_->print(outStream);
  }
};

template<template<typename> class V, typename Real, typename P = Ptr<Vector<Real>>>
inline typename std::enable_if<std::is_base_of<Vector<Real>,V<Real>>::value,P>::type
make_Vector_SimOpt( const Ptr<V<Real>>& vsim, const Ptr<V<Real>>& vopt ) {
  return makePtr<Vector_SimOpt<Real>>(vsim,vopt);
}

} // namespace ROL

#endif
