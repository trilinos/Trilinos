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

#ifndef KOKKOS_CALCULATE_OFFSET_HPP
#define KOKKOS_CALCULATE_OFFSET_HPP

#include <impl/Kokkos_ViewSupport.hpp>
#include <stdint.h>

namespace Kokkos { namespace Impl {

template <typename Layout, typename Shape, typename Enable = void>
struct CalculateOffset;


// LayoutLeft
template <typename Shape>
struct CalculateOffset<LayoutLeft, Shape>
{
  typedef int64_t size_type;
  typedef Shape shape_type;
  typedef LayoutStride<shape_type, LayoutLeft> stride_type;

  //rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply(I0 const& i0, I1 const& i1, stride_type const & stride)
  {
    return i0 + stride.value * i1;
  }

  //rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply(I0 const& i0, I1 const& i1, I2 const& i2, shape_type const& shape, stride_type const & stride)
  {
    return i0 + stride.value * (
           i1 + shape.N1 * i2 );
  }

  //rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i0 + stride.value * (
           i1 + shape.N1 * (
           i2 + shape.N2 * i3 ));
  }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i0 + stride.value * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * i4 )));
  }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i0 + stride.value * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * i5 ))));
  }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i0 + stride.value * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * (
           i5 + shape.N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i0 + stride.value * (
           i1 + shape.N1 * (
           i2 + shape.N2 * (
           i3 + shape.N3 * (
           i4 + shape.N4 * (
           i5 + shape.N5 * (
           i6 + shape.N6 * i7 ))))));
  }

};

// LayoutRight
template <typename Shape>
struct CalculateOffset<LayoutRight, Shape>
{
  typedef int64_t size_type;
  typedef Shape shape_type;
  typedef LayoutStride<shape_type, LayoutRight> stride_type;

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply(I0 const& i0, I1 const& i1, stride_type const & stride)
  {
    return i1 +
           i0 * stride.value;
  }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply(I0 const& i0, I1 const& i1, I2 const& i2, shape_type const& shape, stride_type const & stride)
  {
    return i2 + shape.N2 * ( i1 ) + i0 * stride.value;
  }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i3 + shape.N3 * (
           i2 + shape.N2 * ( i1 )) +
           i0 * stride.value;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * ( i1 ))) +
           i0 * stride.value;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i5 + shape.N5 * (
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * ( i1 )))) +
           i0 * stride.value;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i6 + shape.N6 * (
           i5 + shape.N5 * (
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * ( i1 ))))) +
           i0 * stride.value;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  static size_type apply( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7
                       ,shape_type const& shape, stride_type const & stride)
  {
    return i7 + shape.N7 * (
           i6 + shape.N6 * (
           i5 + shape.N5 * (
           i4 + shape.N4 * (
           i3 + shape.N3 * (
           i2 + shape.N2 * ( i1 )))))) +
           i0 * stride.value;
  }

};

}} // namespace Kokkos::Impl

#endif //KOKKOS_CALCULATE_OFFSET_HPP
