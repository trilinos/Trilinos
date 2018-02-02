
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

#pragma once

#include "XROL_Core.hpp"

namespace XROL {

namespace details {

template<class Element> class StdVector;

template<class Real>
struct IndexType<StdVector<Real>> {
  using type = typename vector<Real>::size_type;
};

template<class Real>
struct ElementType<StdVector<Real>> {
  using type = Real;
}; 

template<class Real>
struct DualType<StdVector<Real>> {
  using type = StdVector<Real>;
};



template<class ElementT>
class StdVector : public Vector<StdVector<ElementT>> {

  // Type aliasing for brevity and consistent labeling
  using IndexT   = index_t<StdVector>;
  using NormT    = norm_t<StdVector>; 
  using DualT    = dual_t<StdVector>;

private: 

  unique_ptr<vector<ElementT>> vec_;

public: 

  // Wraps an existing std::vector
  StdVector( unique_ptr<vector<ElementT>> vec ) : vec_(move(vec)) {}

  // Creates new vector filled with optional value
  StdVector( IndexT dim, ElementT value=0 ) : 
    vec_(move(make_unique<vector<ElementT>>(dim,value))) {}

  void plus( const StdVector& x ) {
    for(IndexT i=0; i<vec_->size(); ++i) (*vec_)[i] += x[i]; 
  }

  void set( const StdVector& x ) {
    for(IndexT i=0; i<vec_->size(); ++i) (*vec_)[i] = x[i]; 
  }

  NormT dot( const StdVector& x ) const {
    ElementT value{0};
    for( IndexT i=0; i<vec_->size(); ++i) value += (*vec_)[i]*x[i];
    return static_cast<NormT>(value);
  }

  NormT norm() const {
    ElementT value{0};
    for( auto e: *vec_ ) value += e*e;
    return sqrt( static_cast<NormT>(value) );
  }

  unique_ptr<Vector<StdVector>> clone() const {
    return move(make_unique<StdVector>(vec_->size()));     
  }

  void axpy( const ElementT alpha, const StdVector& x ) {
    for(size_t i=0; i<vec_->size(); ++i) (*vec_)[i] += alpha*x[i];
  }   

  void fill( const ElementT alpha ) {
    for( auto &e : *vec_ ) e = alpha;
  }

  void scale( const ElementT alpha ) {
    for( auto &e : *vec_ ) e *= alpha;
  }
 
  unique_ptr<Vector<StdVector>> basis( IndexT i ) const {
    auto bp = make_unique<vector<ElementT>>( vec_->size(), 0 );    
    (*bp)[i] = ElementT{1};
    auto b = make_unique<StdVector>( move(bp) );
    return move(b);
  }

  IndexT dimension() const { return vec_->size(); }

  // Does nothing
  const DualT& dual( ) const { 
    return *this; 
  }

  void print( ostream &os, const string& delimiter=" " ) const {
    for( auto e: *vec_ ) os << e << delimiter;
    os << endl;
  } 

  template<class R>
  NormT reduce( R&& r ) const {
    NormT result{r()};
    for( auto e: *vec_ ) result = r(e,result);
    return result;
  }

  template<class F, class... Vs>
  void applyFunction( const F& f, const Vs&... vs ) {
    for( IndexT i=0; i<vec_->size(); ++i ) { 
      (*vec_)[i] = evaluate(f, make_tuple(vs[i]...));
    }
  }
  
  template<class F, class R, class... Vs>
  NormT applyFunctionAndReduce( const F& f, const R& r, const Vs&... vs ) {
    NormT result{r()};
    for( IndexT i=0; i<vec_->size(); ++i ) {
      result = r(result, evaluate(f,make_tuple((*vec_)[i],vs[i]...)));
    }
    return result;
  }


  // Elemental access is well-defined for this type
  // For use in application code without vtable lookup
  ElementT& operator[]( size_t i ) { return (*vec_)[i]; }
  const ElementT& operator[]( size_t i ) const { return (*vec_)[i]; } 

}; // class StdVector

} // namespace details

template<class Real> using StdVector = details::StdVector<Real>;


} // namespace XROL

