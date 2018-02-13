
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
#include <complex>
#include <type_traits>

namespace XROL {
namespace details {

// Scalar Types
template<class T> struct MagnitudeType { using type = T; };
template<class T> struct MagnitudeType<std::complex<T>> { using type = T; };


// Vector Types
template<class V> struct IndexType { using type = int; };

template<class V> struct ElementType { using type = double; };

// Vectors templated on a single type are by default assumed to be 
// templated on element type
template<class E, template<class> class V>
struct ElementType<V<E>> { using type = E; };


// Vectors templated on two types are assumed to have the 
// index and element types in the first and second arguments, respectively
template<class IndexT, class ElementT, template<class,class> class V>
struct IndexType<V<IndexT,ElementT>> { 
  static_assert( std::is_integral<IndexT>::value, "Vectors templated on two parameters "
                 " must have the form Vector<IndexType,ElementType>");
  using type = IndexT; 
};

template<class IndexT, class ElementT, template<class,class> class V>
struct ElementType<V<IndexT,ElementT>> { using type = ElementT; };

template<class V> 
struct NormType { 
  using E    = typename ElementType<V>::type;
  using type = typename MagnitudeType<E>::type; 
};


template<class V> struct DualType { using type = V; };

template<class V> struct IsReflexive {
  static constexpr auto value = is_same_t<V,typename DualType<V>::type>;
};


// Vector Implementation traits
template<class V> 
struct ImplementsElementwise : std::false_type {};


} // namespace details

template<class T> using magnitude_t = typename details::MagnitudeType<T>::type;
template<class V> using index_t     = typename details::IndexType<V>::type;
template<class V> using element_t   = typename details::ElementType<V>::type;
template<class V> using norm_t      = typename details::NormType<V>::type;
template<class V> using dual_t      = typename details::DualType<V>::type;



} // namespace XROL

