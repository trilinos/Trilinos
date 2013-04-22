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

#ifndef KOKKOSARRAY_UTILITY_HPP
#define KOKKOSARRAY_UTILITY_HPP

#include <KokkosArray_Macros.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename T , typename TS1 , typename TS2 = TS1 ,
          bool Enable1 = sizeof(T) == sizeof(TS1) ,
          bool Enable2 = sizeof(T) == sizeof(TS2) >
union UnionPair ;

template< typename T , typename TS2 , bool Enable2 >
union UnionPair<T,T,TS2,true,Enable2> 
{
  typedef T  first_type ;
  typedef T  second_type ;

  first_type  first ;
  second_type second ;

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair() {}

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair( const first_type & rhs ) : first(rhs) {}

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type * cast( first_type * const ptr ) { return ptr ; }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type * cast( const first_type * const ptr ) { return ptr ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type * const ptr ) { return ptr ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type & cast( first_type & ptr ) { return ptr ; }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type & cast( const first_type & ptr ) { return ptr ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type & ptr ) { return ptr ; }
};

template< typename T , typename TS1 , typename TS2 , bool Enable2 >
union UnionPair<T,TS1,TS2,true,Enable2> 
{
  typedef T    first_type ;
  typedef TS1  second_type ;

  first_type  first ;
  second_type second ;

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair() {}

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair( const first_type & rhs ) : first(rhs) {}

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair( const second_type & rhs ) : second(rhs) {}

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type * cast( first_type * const ptr )
  { return reinterpret_cast<second_type*>( ptr ); }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type * cast( const first_type * const ptr )
  { return reinterpret_cast<const second_type*>( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type * const ptr )
  { return reinterpret_cast<volatile second_type*>( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type & cast( first_type & ptr )
  { return reinterpret_cast<second_type&>( ptr ); }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type & cast( const first_type & ptr )
  { return reinterpret_cast<const second_type &>( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type & ptr )
  { return reinterpret_cast<volatile second_type &>( ptr ); }
};

template< typename T , typename TS1 >
union UnionPair<T,TS1,T,false,true> 
{
  typedef T  first_type ;
  typedef T  second_type ;

  first_type  first ;
  second_type second ;

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair() {}

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair( const first_type & rhs ) : first(rhs) {}

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type * cast( first_type * const ptr ) { return ptr ; }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type * cast( const first_type * const ptr ) { return ptr ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type * const ptr ) { return ptr ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type & cast( first_type & ptr ) { return ptr ; }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type & cast( const first_type & ptr ) { return ptr ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type & ptr ) { return ptr ; }
};

template< typename T , typename TS1 , typename TS2 >
union UnionPair<T,TS1,TS2,false,true>
{
  typedef T    first_type ;
  typedef TS2  second_type ;

  first_type  first ;
  second_type second ;

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair() {}

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair( const first_type & rhs ) : first(rhs) {}

  KOKKOSARRAY_INLINE_FUNCTION
  UnionPair( const second_type & rhs ) : second(rhs) {}

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type * cast( first_type * const ptr )
  { return reinterpret_cast<second_type*>( ptr ); }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type * cast( const first_type * const ptr )
  { return reinterpret_cast<const second_type*>( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type * const ptr )
  { return reinterpret_cast<volatile second_type*>( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  second_type & cast( first_type & ptr )
  { return reinterpret_cast<second_type&>( ptr ); }
  
  KOKKOSARRAY_INLINE_FUNCTION
  static
  const second_type & cast( const first_type & ptr )
  { return reinterpret_cast<const second_type &>( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static
  volatile second_type * cast( volatile first_type & ptr )
  { return reinterpret_cast<volatile second_type &>( ptr ); }
};

}
}

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_UTILITY_HPP */

