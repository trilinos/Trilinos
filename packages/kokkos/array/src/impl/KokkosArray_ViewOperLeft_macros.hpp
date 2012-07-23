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

#error "Including <KokkosArray_ViewOperLeft_macros.hpp> without macros defined"

#else

namespace KokkosArray {
namespace Impl {

template< class DataType >
class ViewOperator< DataType , LayoutLeft , KOKKOS_MACRO_DEVICE::memory_space >
{
public:

  typedef typename remove_all_extents< DataType >::type        value_type ;
  typedef typename DefineShape< LayoutLeft , DataType >::type  shape_type ;

private:

  shape_type   m_shape ;
  value_type * m_ptr_on_device ;

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ViewOperator() : m_shape() , m_ptr_on_device(0) {}

  ViewOperator( const ViewOperator & );
  ViewOperator & operator = ( const ViewOperator & );

  template< class , class , class > friend class KokkosArray::View ;
  template< class Dst , class Src >  friend class KokkosArray::Impl::Factory ;

public:

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~ViewOperator() {}

  explicit ViewOperator( const shape_type & shape )
    : m_shape( shape ), m_ptr_on_device(0) {}

  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * ptr_on_device() const { return m_ptr_on_device ; }

  template< typename iType0 >
  inline
  value_type * ptr_on_device( const iType0 & i0 ) const
    {
      typedef typename assert_shape_is_rank_one< shape_type >::type ok_rank ;
      KOKKOS_MACRO_CHECK( Impl::assert_shape_bounds( m_shape, i0 ) );
      return m_ptr_on_device + i0 ;
    }

  template< typename iType0 , typename iType1 >
  inline
  value_type * ptr_on_device( const iType0 & i0 , const iType1 & i1 ) const
    {
      typedef typename assert_shape_is_rank_two< shape_type >::type ok_rank ;
      KOKKOS_MACRO_CHECK( assert_shape_bounds( m_shape, i0, i1 ) );
      return m_ptr_on_device + i0 + m_shape.Stride * i1 ;
    }

  /*------------------------------------------------------------------*/

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator * () const
    {
      typedef typename assert_shape_is_rank_zero< shape_type >::type ok_rank ;
      return *m_ptr_on_device ;
    }

  template< typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator[]( const iType0 & i0 ) const
    {
      typedef typename assert_shape_is_rank_one< shape_type >::type ok_rank ;
      KOKKOS_MACRO_CHECK( assert_shape_bounds( m_shape, i0 ) );
      return m_ptr_on_device[ i0 ];
    }

  /*------------------------------------------------------------------*/

  template< typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ) const
    {
      typedef typename assert_shape_is_rank_one< shape_type >::type ok_rank ;
      KOKKOS_MACRO_CHECK( Impl::assert_shape_bounds( m_shape, i0 ) );
      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      typedef typename assert_shape_is_rank_two< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK( assert_shape_bounds( m_shape, i0, i1 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * i1 ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 ) const
    {
      typedef typename assert_shape_is_rank_three< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * i2 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const
    {
      typedef typename assert_shape_is_rank_four< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * i3 )) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const
    {
      typedef typename assert_shape_is_rank_five< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * i4 ))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const
    {
      typedef typename assert_shape_is_rank_six< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * i5 )))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const
    {
      typedef typename assert_shape_is_rank_seven< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5, i6 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * i6 ))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iType7 & i7 ) const
    {
      typedef typename assert_shape_is_rank_eight< shape_type >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        assert_shape_bounds( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 ) );

      return m_ptr_on_device[ i0 + m_shape.Stride * (
                              i1 + m_shape.N1 * (
                              i2 + m_shape.N2 * (
                              i3 + m_shape.N3 * (
                              i4 + m_shape.N4 * (
                              i5 + m_shape.N5 * (
                              i6 + m_shape.N6 * i7 )))))) ];
    }

#endif /* #if defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

