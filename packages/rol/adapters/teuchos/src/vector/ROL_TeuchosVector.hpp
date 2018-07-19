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

#ifndef ROL_TEUCHOSVECTOR_H
#define ROL_TEUCHOSVECTOR_H

#include <type_traits>

#include "Teuchos_SerialDenseVector.hpp"
#include "ROL_ElementwiseVector.hpp"


/** \class ROL::TeuchosVector 
    \brief Implements the ROL::Vector interface for a Teuchos::SerialDenseVector
 */

namespace ROL {

template<class Ordinal, class Real>
class TeuchosVector : public ElementwiseVector<Real> {

  using SerialDenseVector = Teuchos::SerialDenseVector<Ordinal,Real>;
 
  using UnaryFunction  = ROL::Elementwise::UnaryFunction<Real>;
  using BinaryFunction = ROL::Elementwise::BinaryFunction<Real>;
  using ReductionOp    = ROL::Elementwise::ReductionOp<Real>;

private:

  Teuchos::RCP<SerialDenseVector> vec_;

public:

  TeuchosVector( const Teuchos::RCP<SerialDenseVector> &vec ) : 
    vec_(vec) { }

  TeuchosVector( Ordinal length, bool zeroOut=true ) : 
    vec_( Teuchos::rcp( new SerialDenseVector(length,zeroOut) ) ) {} 

  Real dot( const Vector<Real>& x ) const override { 
    return dynamic_cast<const TeuchosVector&>(x).getVector()->dot(*vec_);
  }

  void plus( const Vector<Real>& x ) override { 
    *vec_ += *(dynamic_cast<const TeuchosVector&>(x).getVector());
  }

  void scale( const Real alpha ) override {
    vec_->scale(alpha);
  }

  int dimension() const override { return static_cast<int>(vec_->length());  }

  void applyUnary( const UnaryFunction &f ) override {
    for( Ordinal i=0; i<vec_->length(); ++i ) {
      (*vec_)(i) = f.apply(((*vec_)(i)));     
    }
  }

  void applyBinary( const BinaryFunction &f, const Vector<Real> &x ) override {
    const TeuchosVector &ex = dynamic_cast<const TeuchosVector&>(x);
    for( Ordinal i=0; i<vec_->length(); ++i ) {
      (*vec_)(i) = f.apply((*vec_)(i),ex(i));
    }    
  }
   
  Real reduce( const ReductionOp &r ) const override {
    Real result = r.initialValue();
    for( Ordinal i=0; i<vec_->length(); ++i ) {
      r.reduce((*vec_)(i),result);
    }
    return result;
  }

  void setScalar(const Real C) override {
    for( Ordinal i=0; i<vec_->length(); ++i ) {
      (*vec_)(i) = C;
    }    
  }

  void randomize(const Real l=0.0, const Real u=1.0) override {
    Real a = (u-l);
    Real b = l;
    Real x(0);
    for( Ordinal i=0; i<vec_->length(); ++i ) {
      x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      (*vec_)(i) = a*x + b;
    }    
  }

  Teuchos::RCP<Vector<Real>> clone() const override { 
    return Teuchos::rcp( new TeuchosVector(vec_->length()) );
  }

  Teuchos::RCP<Vector<Real>> basis( int i ) const override {
    //auto b = Teuchos::rcp( new SerialDenseVector(vec_->length(),true) );
    auto b = Teuchos::rcp( new TeuchosVector(vec_->length(),true) );
    (*b)[static_cast<Ordinal>(i)] = Real{1.0};
    return b;
  }

  Teuchos::RCP<SerialDenseVector> getVector() {
    return vec_;
  }
  
  Teuchos::RCP<const SerialDenseVector> getVector() const { 
    return vec_;
  }
  
  void print( std::ostream &outStream ) const override {
    vec_->print(outStream);
  }

  // Add element access operators to circumvent needing to create and get SerialDenseVectors
  Real& operator() ( Ordinal i ) {
    return (*vec_)(i);
  }

  const Real& operator() ( Ordinal i ) const  {
    return (*vec_)(i); 
  }

  Real& operator[] ( Ordinal i ) {
    return (*vec_)[i];
  }
 
  const Real& operator[] ( Ordinal i ) const {
    return (*vec_)[i]; 
  }

}; // class TeuchosVector


} // namespace ROL


#endif
