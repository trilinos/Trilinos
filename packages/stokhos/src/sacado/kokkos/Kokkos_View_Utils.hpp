// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_UTILS_HPP
#define KOKKOS_VIEW_UTILS_HPP

namespace Kokkos {

namespace Impl {

KOKKOS_INLINE_FUNCTION
void raise_error(const char *msg)
{
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  cuda_abort(msg);
#else
  throw std::runtime_error(msg);
#endif
}

template< class T , class Device > struct RebindStokhosStorageDevice ;

template< class T , class Device >
struct RebindStokhosStorageDevice< T * , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , class Device >
struct RebindStokhosStorageDevice< T [] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , unsigned N , class Device >
struct RebindStokhosStorageDevice< T[N] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type type[N] ;
};

// Get Sacado size from a list of dimensions
template <unsigned Rank> struct GetSacadoSize {};
template <> struct GetSacadoSize<0> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 = 0 ,
                      const size_t n2 = 0 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n0;
  }
};
template <> struct GetSacadoSize<1> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 = 0 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n1;
  }
};
template <> struct GetSacadoSize<2> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n2;
  }
};
template <> struct GetSacadoSize<3> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n3;
  }
};
template <> struct GetSacadoSize<4> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n4;
  }
};
template <> struct GetSacadoSize<5> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n5;
  }
};
template <> struct GetSacadoSize<6> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 ,
                      const size_t n7 = 0 ) {
    return n6;
  }
};
template <> struct GetSacadoSize<7> {
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 ,
                      const size_t n7 ) {
    return n7;
  }
};

template< typename T , T v , bool NonZero = ( v != T(0) ) >
struct integral_nonzero
{
  // Declaration of 'static const' causes an unresolved linker symbol in debug
  // static const T value = v ;
  enum { value = T(v) };
  typedef T value_type ;
  typedef integral_nonzero<T,v> type ;
  KOKKOS_INLINE_FUNCTION integral_nonzero() {}
  KOKKOS_INLINE_FUNCTION integral_nonzero( const T & ) {}
  KOKKOS_INLINE_FUNCTION integral_nonzero( const integral_nonzero & ) {}
  KOKKOS_INLINE_FUNCTION integral_nonzero& operator=(const integral_nonzero &) {return *this;}
};

template< typename T , T zero >
struct integral_nonzero<T,zero,false>
{
  T value ;
  typedef T value_type ;
  typedef integral_nonzero<T,0> type ;
  KOKKOS_INLINE_FUNCTION integral_nonzero() : value() {}
  KOKKOS_INLINE_FUNCTION integral_nonzero( const T & v ) : value(v) {}
  KOKKOS_INLINE_FUNCTION integral_nonzero( const integral_nonzero & v) : value(v.value) {}
  KOKKOS_INLINE_FUNCTION integral_nonzero& operator=(const integral_nonzero & v) { value = v.value; return *this; }
};

} // namespace Impl

} // namespace Kokkos

#endif // KOKKOS_VIEW_UTILS_HPP
