
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

using namespace std;

template<class V>
class Vector : public CRTP<V> {
public:

  using IndexT   = index_t<V>;
  using ElementT = element_t<V>;
  using NormT    = norm_t<V>;
  using DualT    = dual_t<V>;
  
  void plus( const Vector& x )                       { impl().plus(x.impl());        }
  void set( const Vector& x )                        { impl().set(x.impl());         }
  NormT dot( const Vector& x ) const                 { return impl().dot(x.impl());  }
  NormT norm() const                                 { return impl().norm();         }
  unique_ptr<Vector> clone() const                   { return impl().clone();        }
  void axpy( const ElementT alpha, const Vector& x ) { impl().axpy(alpha,x.impl());  }
  void fill( const ElementT alpha )                  { impl().fill(alpha);           }
  void scale( const ElementT alpha )                 { impl().scale(alpha);          }
  unique_ptr basis( IndexT i ) const                 { return impl().basis(i);       }
  IndexT dimension() const                           { return impl().dimension();    }
  void dual(DualT& x) const                          { return impl().dual(x.dual()); }

  void print( ostream& os, const string& delimiter=" " ) const { 
    impl().print(os);            
  } 

  // Elementwise functions

  // y = f(x1,x2,...)
  template<class F, class Vs...>
  void applyFunction( F&& f, Vs&&... vs ) {
    impl().applyFunction( forward<F>(f), forward<Vs>(vs)... );
  }
  
  // result = r(result,y_i) for all i
  template<class R>
  NormT reduce( R&& r ) const {
    impl().reduce(forward<R>(r));
  }

  // result = r(result, f(y_i,x1_i,x2_i,...) for all i
  template<class F, class R, class Vs...>
  NormT applyFunctionAndReduce( F&& f, R&& r, Vs&&... vs ) const {
    return impl().applyFunctionAndReduce( forward<F>(f), forward<R>(r), forward<Vs>(vs)... );  
  }

  


};

template<class V> Vector = details::Vector<V>;

} // namespace details

//template<class V> 
//struct VectorFactory {
//  static std::unique_ptr<V> create() {}
//}; 


} // namespace XROL
