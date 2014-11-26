/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_FUNCTORADAPTER_HPP
#define KOKKOS_FUNCTORADAPTER_HPP

#include <cstddef>
#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Enable = void >
struct FunctorApply
{
  template< class F >
  KOKKOS_FORCEINLINE_FUNCTION static
  void apply( F & ) {}

  template< class F , typename R >
  KOKKOS_FORCEINLINE_FUNCTION static
  void apply( F & , R & ) {}
};

template< class FunctorType >
struct FunctorApply< FunctorType
                   , typename Impl::enable_if< 0 < sizeof( & FunctorType::apply ) >::type >
{
  // use 'enable_if' because FunctorType::apply may or may not have a result argument

  template< class F >
  KOKKOS_FORCEINLINE_FUNCTION static
  typename Impl::enable_if< Impl::is_same< F , FunctorType >::value >::type
    apply( F & f )
      { f.apply(); }

  template< class F , typename R >
  KOKKOS_FORCEINLINE_FUNCTION static
  typename Impl::enable_if< Impl::is_same< F , FunctorType >::value >::type
    apply( F & f , R & r )
      { f.apply(r); }
};

//----------------------------------------------------------------------------

template< class FunctorType , class Enable = void >
struct FunctorInit
{
  template< typename R >
  KOKKOS_FORCEINLINE_FUNCTION static
  void init( const FunctorType & , R & r ) { new(&r) R(); }
};

template< class FunctorType >
struct FunctorInit< FunctorType
                   , typename Impl::enable_if< 0 < sizeof( & FunctorType::init ) >::type >
{
  template< typename R >
  KOKKOS_FORCEINLINE_FUNCTION static
  void init( const FunctorType & f , R & r ) { f.init(r); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_FUNCTORADAPTER_HPP */

