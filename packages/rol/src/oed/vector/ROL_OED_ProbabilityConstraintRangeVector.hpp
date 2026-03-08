// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_PROBABILITYCONSTRAINTRANGEVECTOR_HPP
#define ROL_OED_PROBABILITYCONSTRAINTRANGEVECTOR_HPP

#include "ROL_StdVector.hpp"
#include "ROL_BatchManager.hpp"

namespace ROL::OED {

template<typename Real> 
class ProbabilityConstraintRangeVector : public Vector<Real> {
private:
  const Ptr<StdVector<Real>> vec_;
  const Ptr<BatchManager<Real>> bman_;
  const int dim_;

  const StdVector<Real> & cast(const Vector<Real> &x) const {
    return *static_cast<const ProbabilityConstraintRangeVector<Real>&>(x).getVector();
  }


public:
  ProbabilityConstraintRangeVector(const Ptr<BatchManager<Real>> &bman, int dim, Real val = Real(0))
   : vec_((bman->batchID()==0) ? makePtr<StdVector<Real>>(dim,val) : nullPtr),
     bman_(bman), dim_(dim) {}

  void zero() override {
    if (bman_->batchID()==0) vec_->zero();
    bman_->barrier();
  }

  void set( const Vector<Real> &x ) override {
    if (bman_->batchID()==0) vec_->set(cast(x));
    bman_->barrier();
  }

  void plus( const Vector<Real> &x ) override {
    if (bman_->batchID()==0) vec_->plus(cast(x));
    bman_->barrier();
  }

  void axpy( const Real alpha, const Vector<Real> &x ) override {
    if (bman_->batchID()==0) vec_->axpy(alpha,cast(x));
    bman_->barrier();
  }

  void scale( const Real alpha ) override {
    if (bman_->batchID()==0) vec_->scale(alpha);
    bman_->barrier();
  }

  Real dot( const Vector<Real> &x ) const override {
    Real val(0);
    if (bman_->batchID()==0) val = vec_->dot(cast(x));
    bman_->broadcast(&val,1,0);
    return val;
  }

  Real apply( const Vector<Real> &x ) const override {
    Real val(0);
    if (bman_->batchID()==0) val = vec_->apply(cast(x));
    bman_->broadcast(&val,1,0);
    return val;
  }

  Real norm() const override {
    Real val(0);
    if (bman_->batchID()==0) val = vec_->norm();
    bman_->broadcast(&val,1,0);
    return val;
  }

  Ptr<Vector<Real>> clone() const override {
    return makePtr<ProbabilityConstraintRangeVector<Real>>(bman_,dim_);
  }

  const Vector<Real> & dual() const override {
    return *this;
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) override {
    if (bman_->batchID()==0) vec_->randomize(l,u); 
    bman_->barrier();
  }

  int dimension() const override {
    return dim_;
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) override {
    if (bman_->batchID()==0) vec_->applyUnary(f); 
    bman_->barrier();
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) override {
    if (bman_->batchID()==0) vec_->applyBinary(f,cast(x));
    bman_->barrier();
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const override {
    Real val(0);
    if (bman_->batchID()==0) val = vec_->reduce(r);
    bman_->broadcast(&val,1,0);
    return val;
  }

  void setScalar( const Real C ) override {
    if (bman_->batchID()==0) vec_->setScalar(C); 
    bman_->barrier();
  }

  void print( std::ostream &outStream ) const override {
    if (bman_->batchID()==0) vec_->print(outStream);
    bman_->barrier();
  }

  Ptr<Vector<Real>> basis( const int i ) const override {
    auto e = clone();
    e->setScalar(static_cast<Real>(1));
    return e;
  }

  Real getValue() const {
    Real val(0);
    if (bman_->batchID()==0) val = vec_->getVector()->front();
    bman_->broadcast(&val,1,0);
    return val;
  }

  void setValue(Real val) {
    if (bman_->batchID()==0) (*vec_->getVector())[0] = val;
    bman_->barrier();
  }

  Ptr<const StdVector<Real>> getVector() const {
    if (bman_->batchID()==0) return vec_;
    else return nullPtr;
  }

  Ptr<StdVector<Real>> getVector() {
    if (bman_->batchID()==0) return vec_;
    else return nullPtr;
  }

};

}

#endif
