/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <impl/KokkosArray_Layout_macros.hpp> without macros defined"

#else

#include <impl/KokkosArray_Layout.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template < class T >
struct LayoutMap< Shape<T,Left> , KOKKOS_MACRO_DEVICE::memory_space_new >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space_new::size_type size_type ;
  typedef Shape<T,Left> shape_type ;

  template < typename iType0 , typename iType1 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 )
  {
    return i0 + shape.Stride * ( i1 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 )
  {
    return i0 + shape.Stride * (
           i1 + shape.N1 * i2 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 )
  {
    return i0 + shape.Stride * (
           i1 + shape.N1 * (
           i2 + shape.N2 * i3 ) );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 )
  {
    return i0 + shape.Stride * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * i4 ) ) );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 )
  {
    return i0 + shape.Stride * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * i5 ))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 )
  {
    return i0 + shape.Stride * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * (
           i5 + shape.N5 * i6 )))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 , typename iType7 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 , const iType7 & i7 )
  {
    return i0 + shape.Stride * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * (
           i5 + shape.N5 * (
           i6 + shape.N6 * i7 ))))));
  }
};

//----------------------------------------------------------------------------

template < class T >
struct LayoutMap< Shape<T,Right> , KOKKOS_MACRO_DEVICE::memory_space_new >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space_new::size_type size_type ;
  typedef Shape<T,Right> shape_type ;

  template < typename iType0 , typename iType1 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 )
  {
    return i0 * shape.Stride + i1 ;
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 )
  {
    return i0 * shape.Stride +
           i2 + shape.N2 * (
           i1 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 )
  {
    return i0 * shape.Stride +
           i3 + shape.N3 * (
           i2 + shape.N2 * (
           i1 ) );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 )
  {
    return i0 * shape.Stride +
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * (
           i1 )));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 )
  {
    return i0 * shape.Stride +
           i5 + shape.N5 * (
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * (
           i1 ))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 )
  {
    return i0 * shape.Stride +
           i6 + shape.N6 * (
           i5 + shape.N5 * (
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * (
           i1 )))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 , typename iType7 >
  inline static KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const shape_type & shape ,
                    const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 , const iType7 & i7 )
  {
    return i0 * shape.Stride +
           i7 + shape.N7 * (
           i6 + shape.N6 * (
           i5 + shape.N5 * (
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * (
           i1 ))))));
  }
};

//----------------------------------------------------------------------------

} // Impl namespace
} // KokkosArray namespace

#endif

