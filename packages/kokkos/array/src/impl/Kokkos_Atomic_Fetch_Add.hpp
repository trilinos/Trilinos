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

#if defined( KOKKOSARRAY_ATOMIC_HPP ) && ! defined( KOKKOSARRAY_ATOMIC_FETCH_ADD_HPP )
#define KOKKOSARRAY_ATOMIC_FETCH_ADD_HPP

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined( KOKKOS_ATOMICS_USE_CUDA )

// Support for int, unsigned int, unsigned long long int, and float

KOKKOSARRAY_INLINE_FUNCTION
int atomic_fetch_add( volatile int * const dest , const int val )
{ return atomicAdd((int*)dest,val); }

KOKKOSARRAY_INLINE_FUNCTION
unsigned int atomic_fetch_add( volatile unsigned int * const dest , const unsigned int val )
{ return atomicAdd((unsigned int*)dest,val); }

KOKKOSARRAY_INLINE_FUNCTION
unsigned long long int atomic_fetch_add( volatile unsigned long long int * const dest ,
                                         const unsigned long long int val )
{ return atomicAdd((unsigned long long int*)dest,val); }

KOKKOSARRAY_INLINE_FUNCTION
float atomic_fetch_add( volatile float * const dest , const float val )
{ return atomicAdd((float*)dest,val); }

template < typename T >
KOKKOSARRAY_INLINE_FUNCTION
typename Kokkos::Impl::UnionPair<T,int,unsigned long long int>::first_type
atomic_fetch_add( volatile T * const dest , const T val )
{
  typedef Kokkos::Impl::UnionPair<T,int,unsigned long long int> union_type ;
  typedef typename union_type::second_type type ;

  union_type assumed , old , newval ;

  old.first = *dest ;
  do {
    assumed.second = old.second ;
    newval.first = assumed.first + val ;
    old.second = atomicCAS( (type *) union_type::cast( dest ),
                            assumed.second ,
                            newval.second );
  } while ( assumed.second != old.second );

  return old.first ;
}

//----------------------------------------------------------------------------

#elif defined(KOKKOS_ATOMICS_USE_GCC) || defined(KOKKOS_ATOMICS_USE_INTEL)

KOKKOSARRAY_INLINE_FUNCTION
int atomic_fetch_add( volatile int * const dest , const int val )
{ return __sync_fetch_and_add(dest,val); }

KOKKOSARRAY_INLINE_FUNCTION
long int atomic_fetch_add( volatile long int * const dest , const long int val )
{ return __sync_fetch_and_add(dest,val); }

#if defined( KOKKOS_ATOMICS_USE_GCC )

KOKKOSARRAY_INLINE_FUNCTION
unsigned int atomic_fetch_add( volatile unsigned int * const dest , const unsigned int val )
{ return __sync_fetch_and_add(dest,val); }

KOKKOSARRAY_INLINE_FUNCTION
unsigned long int atomic_fetch_add( volatile unsigned long int * const dest , const unsigned long int val )
{ return __sync_fetch_and_add(dest,val); }

#endif

template < typename T >
KOKKOSARRAY_INLINE_FUNCTION
typename Kokkos::Impl::UnionPair<T,int,long>::first_type
atomic_fetch_add( volatile T * const dest , const T val )
{
  typedef Kokkos::Impl::UnionPair<T,int,long> union_type ;

  union_type assumed , old , newval ;

  old.first = *dest ;
  do {
    assumed.second = old.second ;
    newval.first = assumed.first + val ;
    old.second = __sync_val_compare_and_swap( union_type::cast( dest ),
                                              assumed.second ,
                                              newval.second );
  } while ( assumed.second != old.second );

  return old.first ;
}

//----------------------------------------------------------------------------

#elif defined( KOKKOS_ATOMICS_USE_OMP31 )

template< typename T >
T atomic_fetch_add( volatile T * const dest , const T val )
{
  T retval;
#pragma omp critical
  {
    retval = dest[0];
    dest[0] += val;
  }
  return retval;
}

#endif

//----------------------------------------------------------------------------

}

#endif

