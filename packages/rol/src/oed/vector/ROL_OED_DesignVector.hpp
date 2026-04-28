// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DESIGNVECTOR_HPP
#define ROL_DESIGNVECTOR_HPP

#include "ROL_Vector.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {
namespace OED {

template<typename Real> 
class DesignVector : public Vector<Real> {
private:
  const Ptr<Vector<Real>> vec_;
  Real val_;

  mutable bool isDualInitialized_;
  mutable Ptr<Vector<Real>> dual_vec1_;
  mutable Ptr<DesignVector<Real>> dual_vec_;

public:
  
  // Objective risk only
  DesignVector( const Ptr<Vector<Real>> &vec,
                Real val = Real(0) )
    : vec_(vec), val_(val), isDualInitialized_(false) {}

  void set( const Vector<Real> &x ) {
    const DesignVector<Real> &xs = dynamic_cast<const DesignVector<Real>&>(x);
    vec_->set(*(xs.getVector()));
    val_ = xs.getValue();
  }

  void plus( const Vector<Real> &x ) {
    const DesignVector<Real> &xs = dynamic_cast<const DesignVector<Real>&>(x);
    vec_->plus(*(xs.getVector()));
    val_ += xs.getValue();
  }

  void scale( const Real alpha ) {
    vec_->scale(alpha);
    val_ *= alpha;
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const DesignVector<Real> &xs = dynamic_cast<const DesignVector<Real>&>(x);
    vec_->axpy(alpha,*(xs.getVector()));
    val_ += alpha * xs.getValue();
  }

  Real dot( const Vector<Real> &x ) const {
    const DesignVector<Real> &xs = dynamic_cast<const DesignVector<Real>&>(x);
    return vec_->dot(*(xs.getVector())) + val_*xs.getValue();
  }

  Real norm() const {
    return sqrt( dot(*this) );
  }

  Ptr<Vector<Real>> clone() const {
    return makePtr<DesignVector>(vec_->clone());
  }

  const Vector<Real> &dual() const {
    // Initialize dual vectors if not already initialized
    if ( !isDualInitialized_ ) {
      dual_vec1_ = vec_->dual().clone();
      dual_vec_  = makePtr<DesignVector<Real>>(dual_vec1_);
      isDualInitialized_ = true;
    }
    // Set vector component
    dual_vec1_->set(vec_->dual());
    dual_vec_->setValue(val_);
    // Return dual vector
    return *dual_vec_;
  }

  Ptr<Vector<Real>> basis( const int i )  const {
    Ptr<Vector<Real>> e1;
    Real e2(0);
    int size = vec_->dimension();
    if (i < size && i >= 0) {
      e1 = vec_->basis(i);
    }
    else if (i == size) {
      e1 = vec_->clone();
      e1->zero();
      e2 = 1;
    }
    else {
      ROL_TEST_FOR_EXCEPTION( i < 0 || i > size,
                              std::invalid_argument,
                              "Error: Basis index must be between 0 and vector dimension." );
    }
    return makePtr<DesignVector<Real>>(e1,e2);
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    vec_->applyUnary(f);
    val_ = f.apply(val_);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const DesignVector<Real> &xs = dynamic_cast<const DesignVector<Real>&>(x);
    vec_->applyBinary(f,*xs.getVector());
    val_ = f.apply(val_,xs.getValue());
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    r.reduce(vec_->reduce(r),result);
    r.reduce(val_,result);
    return result;
  }

  void setScalar( const Real C ) {
    vec_->setScalar(C);
    val_ = C;
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    vec_->randomize(l,u);
    val_ = (u-l)*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX) + l;
  }

  int dimension(void) const {
    return vec_->dimension() + 1;
  }

  /***************************************************************************/
  /************ ROL VECTOR ACCESSOR FUNCTIONS ********************************/
  /***************************************************************************/
  Ptr<const Vector<Real>> getVector() const { return vec_; }

  Ptr<Vector<Real>> getVector() { return vec_; }

  Real getValue() const { return val_; }
  void setValue(Real val) { val_ = val; }

};

}
}

#endif
