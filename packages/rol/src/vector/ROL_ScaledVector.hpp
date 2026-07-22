// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_SCALEDVECTOR_HPP
#define ROL_SCALEDVECTOR_HPP

/** @ingroup la_group 
    \class ROL::PrimalScaledVector
    \brief Provides the implementation of the ROL::Vector interface
           that handles scalings in the inner product. A more generic version
           of ROL::PrimalScaledStdVector
*/

/** @ingroup la_group 
    \class ROL::DualScaledVector
    \brief Provides the implementation of the ROL::Vector interface
           that handles scalings in the inner product. A more generic version
           of ROL::PrimalScaledStdVector
*/

#include "ROL_WrappedVector.hpp"
#include "ROL_VectorWorkspace.hpp"

namespace ROL {

// Forward declaration
template<typename Real> class PrimalScaledVector;
template<typename Real> class DualScaledVector;

template<typename Real>
class PrimalScaledVector : public WrappedVector<Real> {

  using V     = Vector<Real>;
  using VPrim = PrimalScaledVector<Real>;
  using VDual = DualScaledVector<Real>;

private:

  mutable Ptv<V>  scaling_vec_;
  mutable VectorWorkspace<Real> workspace_;
  
  Elementwise::Multiply<Real> mult_;

protected:

  VectorWorkspace<Real>& getWorkspace() const { return workspace_; }

  //  y <- y*x elementwise
  void multiply_scaling( const Ptr<V>& y ) const {
    y->applyBinary( mult_, *scaling_vec_ );
  }

public: 

  PrimalScaledVector( const Ptr<V>& vec, const Ptr<V>& scaling_vec ) :
    WrappedVector<Real>(vec), scaling_vec_(scaling_vec) {}

  virtual ~PrimalScaledVector() {}

  virtual Real dot( const V& x ) const override {
    auto y = workspace_.copy(x);
    multiply_scaling( y );  
    return this->getVector()->dot(*y);
  }

  virtual Ptr<V> clone() const override { 
    return makePtr<VPrim>( this->getVector()->clone(), scaling_vec_ );
  }

  virtual Ptr<V> basis( const int i ) const override {
    return makePtr<VPrim>( this->getVector()->basis(i), scaling_vec_ );
  }

  virtual void const V& dual() const override {
    auto dual_vec = workspace_.copy( this->getVector() );
    multiply_scaling( dual_vec );
    return *( makePtr<VDual>( dual_vec, scaling_vec ) );
  } 

  const Ptr<V>& getScalingVector() { return scaling_vec_; }
  const Ptr<const V>& getScalingVector() const { return scaling_vec_; }

  void setScalingVector( const Ptr<const V&>& scaling_vec ) const { 
    scaling_vec_ = scaling_vec;
  }

}; // class PrimalScaledVector



template<typename Real>
class DualScaledVector : public WrappedVector<Real> {

  using V     = Vector<Real>;
  using VPrim = PrimalScaledVector<Real>;
  using VDual = DualScaledVector<Real>;

private:

  mutable Ptv<V>  scaling_vec_;
  mutable VectorWorkspace<Real> workspace_;
  
  Elementwise::Divide<Real> div_;

protected:

  VectorWorkspace<Real>& getWorkspace() const { return workspace_; }

  //  y <- y/x elementwise
  void divide_scaling( const <V>& y ) const {
    y->applyBinary( div_, *scaling_vec_ );
  }

public: 

  DualScaledVector( const Ptr<V>& vec, const Ptr<V>& scaling_vec ) :
    WrappedVector<Real>(vec), scaling_vec_(scaling_vec) {}

  virtual ~DualScaledVector() {}

  virtual Real dot( const V& x ) const override {
    auto y = workspace_.copy(x);
    divide_scaling( y );
    return this->getVector()->dot(*y);
  }

  virtual Ptr<V> clone() const override { 
    return makePtr<VDual>( this->getVector()->clone(), scaling_vec_ );
  }

  virtual Ptr<V> basis( const int i ) const override {
    return makePtr<VDual>( this->getVector()->basis(i), scaling_vec_ );
  }

  virtual void const V& dual() const override {
    auto primal_vec = workspace_.copy( this->getVector() );
    divide_scaling( primal_vec );
    return *( makePtr<VPrim>( primal_vec, scaling_vec ) );
  } 

  const Ptr<V>& getScalingVector() { return scaling_vec_; }
  const Ptr<const V>& getScalingVector() const { return scaling_vec_; }

  void setScalingVector( const Ptr<const V&>& scaling_vec ) const { 
    scaling_vec_ = scaling_vec;
  }

}; // class PrimalScaledVector

} // namespace ROL


#endif
