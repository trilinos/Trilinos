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

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <KokkosArray_ViewOperRight_macros.hpp> without macros defined"

#else

#include <impl/KokkosArray_ShapeRight.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , 0 , 0 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                      * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,0,0> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator * () const
    { return *m_ptr_on_device ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( const ValueType & v ) const
    { *m_ptr_on_device = v ; }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 1 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,1> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator[]( const iType0 & i0 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 2 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,2> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1 ) );

      return m_ptr_on_device[ i1 + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 3 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,3> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2 ) );

      return m_ptr_on_device[ i2 + m_shape.N2 * (
                              i1 ) + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 4 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,4> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3 ) );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 5 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,5> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4 ) );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 6 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,6> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5 ) );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 7 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,7> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ,
                          const iType6 & i6 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5, i6 ) );

      return m_ptr_on_device[ i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))))) + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

template< typename ValueType , class ShapeType , unsigned RankDyn >
class ViewOper< ValueType , KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                Shape< LayoutRight , ShapeType , RankDyn , 8 > >
{
private:
  template< class , class , class > friend class KokkosArray::View ;

  ValueType                            * m_ptr_on_device ;
  Shape<LayoutRight,ShapeType,RankDyn,8> m_shape ;

public:

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  ValueType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ,
                          const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOSARRAY_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 ) );

      return m_ptr_on_device[ i7 + m_shape.N7 * (
                              i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))))) + i0 * m_shape.Stride ];
    }

#endif

};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

