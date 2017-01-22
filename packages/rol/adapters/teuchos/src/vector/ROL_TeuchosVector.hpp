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

#include "Teuchos_SerialDenseVector.hpp"
#include "ROL_ElementwiseVector.hpp"


/** \class ROL::TeuchosVector 
    \brief Implements the ROL::Vector interface for a Teuchos::SerialDenseVector
 */

namespace ROL {

template<class Ordinal, class Real>
class TeuchosVector : public ElementwiseVector<Real> {

  template <typename T> using RCP = Teuchos::RCP<T>;

  typedef Teuchos::SerialDenseVector<Ordinal,Real>  SDV;
  typedef TeuchosVector<Ordinal,Real>               TV;
  typedef Vector<Real>                              V;

  typedef ROL::Elementwise::UnaryFunction<Real>     UF;
  typedef ROL::Elementwise::BinaryFunction<Real>    BF;
  typedef ROL::Elementwise::ReductionOp<Real>       RO;

private:

  RCP<SDV> vec_;
  Ordinal dim_;

public:

  TeuchosVector( const RCP<SDV> &vec ) : vec_(vec), dim_(vec_->length()) { }

  // Create a vector of given length with optional zeroing
  TeuchosVector( Ordinal length, bool zeroOut=true ) :
    vec_(Teuchos::rcp( new SDV(length,zeroOut) ) ),
    dim_(length) {
  }
  
  int dimension() {
    return static_cast<int>(dim_);
  }

  RCP<V> basis( const int i ) const {
    RCP<TV> b = Teuchos::rcp( new TV(dim_,true) );
    (*b)[static_cast<Ordinal>(i)] = Real(1.0);
    return b;
  }

  void applyUnary( const UF &f ) {
  
    for( Ordinal i=0; i<dim_; ++i ) {
      (*vec_)(i) = f.apply(((*vec_)(i)));     
    }
  }

  void applyBinary( const BF &f, const V &x ) {
  
//    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
//                                std::invalid_argument,
//                                "Error: Vectors must have the same dimension." );

    const TV &ex = Teuchos::dyn_cast<const TV>(x);
    for( Ordinal i=0; i<dim_; ++i ) {
      (*vec_)(i) = f.apply((*vec_)(i),ex(i));
    }    
  }
   
  Real reduce( const RO &r ) const {
    Real result = r.initialValue();
    for( Ordinal i=0; i<dim_; ++i ) {
      r.reduce((*vec_)(i),result);
    }
    return result;
  }

  RCP<V> clone() const {
    return Teuchos::rcp( new TV( dim_ ) );
  }

  RCP<SDV> getVector() {
    return vec_;
  }
  
  RCP<const SDV> getVector() const { 
    return vec_;
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
