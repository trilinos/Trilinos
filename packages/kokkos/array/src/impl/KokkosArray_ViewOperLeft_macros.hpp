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

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <KokkosArray_ViewOperLeft_macros.hpp> without macros defined"

#else

#include <impl/KokkosArray_ShapeLeft.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType , Shape< LayoutLeft , ValueSize, 0 > >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType                   * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,0> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator * () const
    { return *m_ptr_on_device ; }

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const ValueType & v ) const
    { *m_ptr_on_device = v ; }

  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( void ) const
    { return *m_ptr_on_device ; }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize , unsigned s0 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType , Shape< LayoutLeft, ValueSize, 1, s0 > >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType                      * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,1,s0> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator[]( const iType0 & i0 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize , unsigned s0 , unsigned s1 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType , Shape< LayoutLeft, ValueSize, 2, s0, s1 > >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType                         * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,2,s0,s1> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * i1 ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize ,
          unsigned s0 , unsigned s1 , unsigned s2 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType ,
                Shape<LayoutLeft,ValueSize,3,s0,s1,s2> >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,3,s0,s1,s2> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * i2 ) ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType ,
                Shape<LayoutLeft,ValueSize,4,s0,s1,s2,s3> >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,4,s0,s1,s2,s3> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType ,
                Shape<LayoutLeft,ValueSize,5,s0,s1,s2,s3,s4> >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,5,s0,s1,s2,s3,s4> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * i4 ))) ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType ,
                Shape<LayoutLeft,ValueSize,6,s0,s1,s2,s3,s4,s5> >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,6,s0,s1,s2,s3,s4,s5> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * i5 )))) ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType ,
                Shape<LayoutLeft,ValueSize,7,s0,s1,s2,s3,s4,s5,s6> >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,7,s0,s1,s2,s3,s4,s5,s6> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ,
                          const iType6 & i6 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5, i6 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * i6 ))))) ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , unsigned ValueSize ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
class ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                ValueType ,
                Shape<LayoutLeft,ValueSize,8,s0,s1,s2,s3,s4,s5,s6,s7> >
{
private:
  template< class , class , class , class > friend class KokkosArray::View ;

  ValueType * m_ptr_on_device ;
  Shape<LayoutLeft,ValueSize,8,s0,s1,s2,s3,s4,s5,s6,s7> m_shape ;

public:

#if defined( KOKKOSARRAY_INLINE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  KOKKOSARRAY_INLINE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ,
                          const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * (
                              i6 + m_shape.N6 * i7 )))))) ];
    }

#endif

};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

