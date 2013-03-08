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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_VIEWASSIGNMENT_HPP
#define KOKKOSARRAY_VIEWASSIGNMENT_HPP

#include <typeinfo>
#include <utility>
#include <KokkosArray_Macros.hpp>

#include <impl/KokkosArray_ArrayTraits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class DstViewType >
struct ViewInitialize { static void apply( const DstViewType & ) {} };

template< class T , class L , class D , class M , class S >
size_t allocation_count( const View<T,L,D,M,S> & view )
{
  return ViewAssignment<S>::allocation_count( view );
}

//----------------------------------------------------------------------------

template< class ViewTraits ,
          class MemorySpace  = typename ViewTraits::memory_space ,
          class MemoryTraits = typename ViewTraits::memory_traits ,
          class ExecSpec     = KokkosArray::ExecutionSpace >
struct ViewTracking {
  KOKKOSARRAY_INLINE_FUNCTION static void increment( const void * ) {}
  KOKKOSARRAY_INLINE_FUNCTION static void decrement( const void * ) {}
};

template< class ViewTraits , class MemorySpace , class MemoryTraits >
struct ViewTracking< ViewTraits , MemorySpace , MemoryTraits ,
          typename enable_if< MemoryTraits::managed , HostSpace >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION static void increment( const void * ptr )
    { MemorySpace::increment( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION static void decrement( const void * ptr )
    { MemorySpace::decrement( ptr ); }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstView , class SrcView ,
          class DstValueType  = typename DstView::value_type ,
          class DstValueSpace = typename DstView::memory_space ,
          class SrcValueType  = typename SrcView::value_type ,
          class SrcValueSpace = typename SrcView::memory_space >
struct ValueCompatible ;

template< class DstView , class SrcView , class ValueType , class ValueSpace >
struct ValueCompatible< DstView , SrcView ,
                        ValueType , ValueSpace ,
                        ValueType , ValueSpace >
{
  typedef ValueType type ;
  enum { value = true };
};

template< class DstView , class SrcView , class ValueType , class ValueSpace >
struct ValueCompatible< DstView , SrcView ,
                        const ValueType , ValueSpace ,
                              ValueType , ValueSpace >
{
  typedef ValueType type ;
  enum { value = true };
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_VIEWASSIGNMENT_HPP */

