/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#if defined( KOKKOSARRAY_ATOMIC_HPP ) && ! defined( KOKKOSARRAY_ATOMIC_COMPARE_EXCHANGE_STRONG_HPP )
#define KOKKOSARRAY_ATOMIC_COMPARE_EXCHANGE_STRONG_HPP

namespace Kokkos {

//----------------------------------------------------------------------------
// Cuda native CAS supports int, unsigned int, and unsigned long long int (non-standard type).
// Must cast-away 'volatile' for the CAS call.

#if defined( KOKKOS_ATOMICS_USE_CUDA )

KOKKOSARRAY_INLINE_FUNCTION
int atomic_compare_exchange( volatile int * const dest, const int compare, const int val)
{ return atomicCAS((int*)dest,compare,val); }

KOKKOSARRAY_INLINE_FUNCTION
unsigned int atomic_compare_exchange( volatile unsigned int * const dest, const unsigned int compare, const unsigned int val)
{ return atomicCAS((unsigned int*)dest,compare,val); }

KOKKOSARRAY_INLINE_FUNCTION
unsigned long long int atomic_compare_exchange( volatile unsigned long long int * const dest ,
                                                const unsigned long long int compare ,
                                                const unsigned long long int val )
{ return atomicCAS((unsigned long long int*)dest,compare,val); }

template < typename T >
KOKKOSARRAY_INLINE_FUNCTION
typename Kokkos::Impl::UnionPair<T,int,unsigned long long int>::first_type
atomic_compare_exchange( volatile T * const dest , const T compare , const T val )
{
  typedef Kokkos::Impl::UnionPair<T,int,unsigned long long int> union_type ;
  typedef typename union_type::second_type int_type ;

  return union_type( atomicCAS( (int_type *) union_type::cast( dest ) ,
                                union_type::cast( compare ) ,
                                union_type::cast( val ) )
                   ).first ;
}

//----------------------------------------------------------------------------
// GCC native CAS supports int, long, unsigned int, unsigned long.
// Intel native CAS support int and long with the same interface as GCC.

#elif defined(KOKKOS_ATOMICS_USE_GCC) || defined(KOKKOS_ATOMICS_USE_INTEL)

KOKKOSARRAY_INLINE_FUNCTION
int atomic_compare_exchange( volatile int * const dest, const int compare, const int val)
{ return __sync_val_compare_and_swap(dest,compare,val); }

KOKKOSARRAY_INLINE_FUNCTION
long atomic_compare_exchange( volatile long * const dest, const long compare, const long val )
{ return __sync_val_compare_and_swap(dest,compare,val); }

#if defined( KOKKOS_ATOMICS_USE_GCC )

// GCC supports unsigned

KOKKOSARRAY_INLINE_FUNCTION
unsigned int atomic_compare_exchange( volatile unsigned int * const dest, const unsigned int compare, const unsigned int val )
{ return __sync_val_compare_and_swap(dest,compare,val); }

KOKKOSARRAY_INLINE_FUNCTION
unsigned long atomic_compare_exchange( volatile unsigned long * const dest ,
                                       const unsigned long compare ,
                                       const unsigned long val )
{ return __sync_val_compare_and_swap(dest,compare,val); }

#endif

template < typename T >
KOKKOSARRAY_INLINE_FUNCTION
typename Kokkos::Impl::UnionPair<T,int,long>::first_type
atomic_compare_exchange( volatile T * const dest, const T compare, const T val )
{
  typedef Kokkos::Impl::UnionPair<T,int,long> union_type ;

  return union_type(
    __sync_val_compare_and_swap( union_type::cast( dest ) ,
                                 union_type::cast( compare ) ,
                                 union_type::cast( val ) )
  ).first ;
}

//----------------------------------------------------------------------------

#elif defined( KOKKOS_ATOMICS_USE_OMP31 )

template< typename T >
KOKKOSARRAY_INLINE_FUNCTION
T atomic_compare_exchange( volatile T * const dest, const T compare, const T val )
{
  T retval;
#pragma omp critical
  {
    retval = dest[0];
    if ( retval == compare )
  	dest[0] = val;
  }
  return retval;
}

#endif

//----------------------------------------------------------------------------

} // namespace Kokkos

#endif

