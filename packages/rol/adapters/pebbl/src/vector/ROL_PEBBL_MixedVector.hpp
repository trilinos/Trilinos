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

#ifndef ROL_PEBBL_MIXEDVECTOR_H
#define ROL_PEBBL_MIXEDVECTOR_H

#include "ROL_Vector.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::MixedVector
    \brief Defines the pebbl mixed vector interface.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Real>
class MixedVector : public Vector<Real> {
private:
  Ptr<Vector<Real>> cntvec_;
  Ptr<Vector<Real>> intvec_;
  mutable ROL::Ptr<Vector<Real>> dual_cntvec_;
  mutable ROL::Ptr<Vector<Real>> dual_intvec_;
  mutable ROL::Ptr<MixedVector<Real>> dual_vec_;

public:
  MixedVector(const Ptr<Vector<Real>> &cntvec,
              const Ptr<Vector<Real>> &intvec)
    : cntvec_(cntvec), intvec_(intvec) {
    dual_cntvec_ = cntvec_->dual().clone();
    dual_intvec_ = intvec_->dual().clone();
  }

  MixedVector(const Ptr<Vector<Real>> &intvec)
    : cntvec_(nullPtr), intvec_(intvec) {
    dual_intvec_ = intvec_->dual().clone();
  }

  void plus( const Vector<Real> &x ) {
    const MixedVector<Real> &xs = dynamic_cast<const MixedVector<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    if (cntvec_ != nullPtr) cntvec_->plus(*(xs.getContinuousVariables()));
    intvec_->plus(*(xs.getIntegerVariables()));
  }   

  void scale( const Real alpha ) {
    if (cntvec_ != nullPtr) cntvec_->scale(alpha);
    intvec_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const MixedVector<Real> &xs = dynamic_cast<const MixedVector<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    if (cntvec_ != nullPtr) cntvec_->axpy(alpha,*(xs.getContinuousVariables()));
    intvec_->axpy(alpha,*(xs.getIntegerVariables()));
  }

  Real dot( const Vector<Real> &x ) const {
    const MixedVector<Real> &xs = dynamic_cast<const MixedVector<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    Real prod(0);
    if (cntvec_ != nullPtr) prod = cntvec_->dot(*(xs.getContinuousVariables()));
    return prod + intvec_->dot(*(xs.getIntegerVariables()));
  }

  Real norm() const {
    Real norm1(0);
    if (cntvec_ != nullPtr) norm1 = cntvec_->norm();
    Real norm2 = intvec_->norm();
    return sqrt( norm1*norm1 + norm2*norm2 );
  } 

  ROL::Ptr<Vector<Real> > clone() const {
    if (cntvec_ != nullPtr) 
      return ROL::makePtr<MixedVector>(cntvec_->clone(),intvec_->clone());  
    return ROL::makePtr<MixedVector>(intvec_->clone());  
  }

  const Vector<Real> & dual(void) const {
    dual_intvec_->set(intvec_->dual());
    if (cntvec_ != nullPtr) {
      dual_cntvec_->set(cntvec_->dual());
      dual_vec_ = ROL::makePtr<MixedVector<Real>>(dual_cntvec_,dual_intvec_); 
    }
    else {
      dual_vec_ = ROL::makePtr<MixedVector<Real>>(dual_intvec_); 
    }
    return *dual_vec_;
  }

  Real apply( const Vector<Real> &x ) const {
    const MixedVector<Real> &xs = dynamic_cast<const MixedVector<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    Real prod(0);
    if (cntvec_ != nullPtr) prod = cntvec_->apply(*(xs.getContinuousVariables()));
    return prod + intvec_->apply(*(xs.getIntegerVariables()));
  }

  ROL::Ptr<Vector<Real> > basis( const int i )  const {
    int n1(0);
    if (cntvec_ != nullPtr) n1 = (cntvec_)->dimension();
    if ( i < n1 ) {
      ROL::Ptr<Vector<Real>> e1 = (cntvec_)->basis(i);
      ROL::Ptr<Vector<Real>> e2 = (intvec_)->clone(); e2->zero();
      ROL::Ptr<Vector<Real>> e  = ROL::makePtr<MixedVector<Real>>(e1,e2);
      return e;
    }
    else {
      ROL::Ptr<Vector<Real>> e1 = (cntvec_)->clone(); e1->zero();
      ROL::Ptr<Vector<Real>> e2 = (intvec_)->basis(i-n1);
      ROL::Ptr<Vector<Real>> e  = ROL::makePtr<MixedVector<Real>>(e1,e2);
      return e;
    }
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    if (cntvec_ != nullPtr) cntvec_->applyUnary(f);
    intvec_->applyUnary(f);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const MixedVector<Real> &xs = dynamic_cast<const MixedVector<Real>&>(x);
    if (cntvec_ != nullPtr) cntvec_->applyBinary(f,*xs.getContinuousVariables());
    intvec_->applyBinary(f,*xs.getIntegerVariables());
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    if (cntvec_ != nullPtr) r.reduce(cntvec_->reduce(r),result);
    r.reduce(intvec_->reduce(r),result);
    return result;
  }

  void setScalar( const Real C ) {
    if (cntvec_ != nullPtr) cntvec_->setScalar(C);
    intvec_->setScalar(C);
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    if (cntvec_ != nullPtr) cntvec_->randomize(l,u);
    intvec_->randomize(l,u);
  }


  int dimension() const {
    int n1(0);
    if (cntvec_ != nullPtr) n1 = cntvec_->dimension();
    return n1 + (intvec_)->dimension();
  }

  ROL::Ptr<const Vector<Real>> getContinuousVariables() const { 
    return cntvec_; 
  }

  ROL::Ptr<const Vector<Real>> getIntegerVariables() const { 
    return intvec_; 
  }

  ROL::Ptr<Vector<Real>> getContinuousVariables() { 
    return cntvec_; 
  }

  ROL::Ptr<Vector<Real>> getIntegerVariables() { 
    return intvec_; 
  }

  void setContinuousVariables(const Vector<Real>& vec) { 
    if (cntvec_ != nullPtr) cntvec_->set(vec);
  }
  
  void setIntegerVariables(const Vector<Real>& vec) { 
    intvec_->set(vec); 
  }

  void print( std::ostream &outStream ) const {
    if (cntvec_ != nullPtr) {
      outStream << "Continuous: ";
      cntvec_->print(outStream);
    }
    outStream << "Integer: ";
    intvec_->print(outStream);
  }
}; // class MixedVector

} // namespace PEBBL
} // namespace ROL

#endif
