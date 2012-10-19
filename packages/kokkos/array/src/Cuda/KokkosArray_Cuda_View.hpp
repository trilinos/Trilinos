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

#include <KokkosArray_Host.hpp>
#include <KokkosArray_View.hpp>

#include <impl/KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_ViewOperLeft_macros.hpp>
#include <impl/KokkosArray_ViewOperRight_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <impl/KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType >
void View< DataType , LayoutType , Cuda >::internal_private_create(
  const std::string & label ,
  const typename View< DataType , LayoutType , Cuda >::shape_type shape )
{
  typedef typename View< DataType , LayoutType , Cuda >::shape_type shape_type ;

  const size_t count = Impl::ShapeMap<shape_type>::allocation_count( shape );

  oper_type::m_shape = shape ;
  oper_type::m_ptr_on_device = (scalar_type *)
    memory_space::allocate( label ,
                            typeid(scalar_type) ,
                            sizeof(scalar_type) ,
                            count );
}

//----------------------------------------------------------------------------

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

