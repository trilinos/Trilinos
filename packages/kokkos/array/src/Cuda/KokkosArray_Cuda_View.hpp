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

#ifndef KOKKOSARRAY_CUDA_VIEW_HPP
#define KOKKOSARRAY_CUDA_VIEW_HPP

#include <KokkosArray_View.hpp>
#include <KokkosArray_HostSpace.hpp>
#include <KokkosArray_CudaSpace.hpp>
#include <Cuda/KokkosArray_Cuda_abort.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template<>
struct AssertShapeBoundsAbort< CudaSpace >
{
  KOKKOSARRAY_INLINE_FUNCTION
  static void apply( const size_t /* rank */ ,

                     const size_t /* n0 */ , const size_t /* n1 */ ,
                     const size_t /* n2 */ , const size_t /* n3 */ ,
                     const size_t /* n4 */ , const size_t /* n5 */ ,
                     const size_t /* n6 */ , const size_t /* n7 */ ,

                     const size_t /* i0 */ , const size_t /* i1 */ ,
                     const size_t /* i2 */ , const size_t /* i3 */ ,
                     const size_t /* i4 */ , const size_t /* i5 */ ,
                     const size_t /* i6 */ , const size_t /* i7 */ )
    {
      KokkosArray::cuda_abort("KokkosArray::View array bounds violation");
    }
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename ValueType , class LayoutSrc >
inline
void deep_copy( ValueType & dst ,
                const View< ValueType , LayoutSrc , Cuda > & src )
{
  typedef View< ValueType , LayoutSrc , Cuda > src_type ;
  typedef typename src_type::shape_type        src_shape ;

  typedef typename Impl::assert_shape_is_rank_zero< src_shape >::type ok_rank ;

  DeepCopy<HostSpace,CudaSpace>( & dst , src.ptr_on_device() , sizeof(ValueType) );
}

template< typename ValueType , class LayoutDst >
inline
void deep_copy( const View< ValueType , LayoutDst , Cuda > & dst ,
                const ValueType & src )
{
  typedef View< ValueType , LayoutDst , Cuda > dst_type ;
  typedef typename dst_type::shape_type        dst_shape ;

  typedef typename Impl::assert_shape_is_rank_zero< dst_shape >::type ok_rank ;

  DeepCopy<CudaSpace,HostSpace>( dst.ptr_on_device() , & src , sizeof(ValueType) );
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_CUDA_VIEW_HPP */

