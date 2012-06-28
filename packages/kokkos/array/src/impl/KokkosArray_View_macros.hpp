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

#error "Including <KokkosArray_View_macros.hpp> without macros defined"

#else

namespace KokkosArray {

template< class DataType , class LayoutSpec >
class View< DataType , LayoutSpec , KOKKOS_MACRO_DEVICE >
{
public:
  typedef View< DataType , LayoutSpec , KOKKOS_MACRO_DEVICE >  type ;
  typedef View< DataType , LayoutSpec , Host >  HostMirror ;

  typedef DataType             data_type ;
  typedef LayoutSpec           layout_spec ;
  typedef KOKKOS_MACRO_DEVICE  device_type ;

  typedef typename Impl::remove_all_extents<data_type>::type  value_type ;
  typedef typename LayoutSpec::array_layout                   array_layout ;
  typedef typename device_type::memory_space_new              memory_space ;
  typedef typename device_type::size_type                     size_type ;

private:

  typedef typename Impl::change_empty_extent_to_zero_extent< 
          typename Impl::remove_const< data_type >::type >::type
    clean_data_type ;

  typedef Impl::Shape< array_layout , clean_data_type > shape_type ;

  typedef Impl::ShapeMap< shape_type , memory_space > shape_map ;

public:

  /*------------------------------------------------------------------*/

  static const unsigned Rank = shape_type::rank ;

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return shape_type::rank ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_0() const { return m_shape.N0 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_1() const { return m_shape.N1 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_2() const { return m_shape.N2 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_3() const { return m_shape.N3 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_4() const { return m_shape.N4 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_5() const { return m_shape.N5 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_6() const { return m_shape.N6 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_7() const { return m_shape.N7 ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & r ) const
    {
      return iType(0) == r ? m_shape.N0 : (
             iType(1) == r ? m_shape.N1 : (
             iType(2) == r ? m_shape.N2 : (
             iType(3) == r ? m_shape.N3 : (
             iType(4) == r ? m_shape.N4 : (
             iType(5) == r ? m_shape.N5 : (
             iType(6) == r ? m_shape.N6 : (
             iType(7) == r ? m_shape.N7 : 0 )))))));
    }

  /** \brief  Because memory is contiguous this is exposed */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * ptr_on_device() const
  { return m_ptr_on_device ; }

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator * () const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 0 >::type ok_rank ;
      return *m_ptr_on_device ;
    }

  template< typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator[]( const iType0 & i0 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 1 >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ]; 
    }

  template< typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 1 >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ]; 
    }

  template< typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 2 >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0, i1 ) );

      return m_ptr_on_device[ shape_map::offset( m_shape, i0, i1 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 3 >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0, i1, i2 ) );

      return m_ptr_on_device[ shape_map::offset( m_shape, i0 , i1 , i2 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 4 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3 ) );

      return m_ptr_on_device[ shape_map::offset( m_shape, i0, i1, i2, i3 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 5 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4 ) );

      return m_ptr_on_device[
               shape_map::offset( m_shape, i0, i1, i2, i3, i4 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 6 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5 ) );

      return m_ptr_on_device[
               shape_map::offset( m_shape, i0, i1, i2, i3, i4, i5 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 7 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5, i6 ) );

      return m_ptr_on_device[
               shape_map::offset( m_shape, i0, i1, i2, i3, i4, i5, i6 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ,
                           const iType7 & i7 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 8 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 ) );

      return m_ptr_on_device[ shape_map::offset(m_shape,i0,i1,i2,i3,i4,i5,i6,i7) ];
    }

#endif /* #if defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  template< typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 1 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0) );
    }

  template< typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 2 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1) );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ,
          const iType2 & i2 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 3 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1,i2) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ,
          const iType2 & i2 , const iType3 & i3 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 4 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1,i2,i3) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ,
          const iType2 & i2 , const iType3 & i3 ,
          const iType4 & i4 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 5 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1,i2,i3,i4) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ,
          const iType2 & i2 , const iType3 & i3 ,
          const iType4 & i4 , const iType5 & i5 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 6 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1,i2,i3,i4,i5) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ,
          const iType2 & i2 , const iType3 & i3 ,
          const iType4 & i4 , const iType5 & i5 ,
          const iType6 & i6 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 7 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5, i6 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1,i2,i3,i4,i5,i6) );
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View<value_type,device_type>
    view( const iType0 & i0 , const iType1 & i1 ,
          const iType2 & i2 , const iType3 & i3 ,
          const iType4 & i4 , const iType5 & i5 ,
          const iType6 & i6 , const iType7 & i7 ) const
    {
      typedef typename Impl::assert_shape_is_rank< shape_type , 8 >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 ) );

      return View<value_type,device_type>(
             m_ptr_on_device + shape_map::offset(m_shape,i0,i1,i2,i3,i4,i5,i6,i7) );
    }

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View()
    : m_shape()
    , m_ptr_on_device(0)
    {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View( const View & rhs )
    : m_shape(         rhs.m_shape )
    , m_ptr_on_device( rhs.m_ptr_on_device )
    { memory_space::increment( m_ptr_on_device ); }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View & operator = ( const View & rhs )
  {
    memory_space::decrement( m_ptr_on_device );
    m_shape          = rhs.m_shape ;
    m_ptr_on_device  = rhs.m_ptr_on_device ;
    memory_space::increment( m_ptr_on_device );
    return *this ;
  }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~View() { memory_space::decrement( m_ptr_on_device ); }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a compatible view */

  template< class rhsType , class rhsMapSpec , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View( const View< rhsType , rhsMapSpec , rhsMemory > & rhs )
    : m_shape(         rhs.m_shape ) // Must be same type
    , m_ptr_on_device( rhs.m_ptr_on_device ) // preserves 'const' requirement
    {
      typedef View< rhsType , rhsMapSpec , rhsMemory > rhs_type ;
      typedef typename rhs_type::array_layout rhs_array_layout ;
      typedef typename rhs_type::memory_space rhs_memory_space ;
      typedef typename rhs_type::value_type   rhs_value_type ;

      typedef typename Impl::StaticAssertAssignable<value_type,  rhs_value_type>  ::type ok_value ;
      typedef typename Impl::StaticAssertSame<array_layout,rhs_array_layout>::type ok_layout ;
      typedef typename Impl::StaticAssertSame<memory_space,rhs_memory_space>::type ok_memory ;

      memory_space::increment( m_ptr_on_device );
    }

  template< class rhsType , class rhsMapSpec , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View & operator = ( const View< rhsType , rhsMapSpec , rhsMemory > & rhs )
    {
      typedef View< rhsType , rhsMapSpec , rhsMemory > rhs_type ;
      typedef typename rhs_type::array_layout rhs_array_layout ;
      typedef typename rhs_type::memory_space rhs_memory_space ;
      typedef typename rhs_type::value_type   rhs_value_type ;

      typedef typename Impl::StaticAssertAssignable<value_type,  rhs_value_type>  ::type ok_value ;
      typedef typename Impl::StaticAssertSame<array_layout,rhs_array_layout>::type ok_layout ;
      typedef typename Impl::StaticAssertSame<memory_space,rhs_memory_space>::type ok_memory ;

      memory_space::decrement( m_ptr_on_device );
      m_shape          = rhs.m_shape ; // Must be same type
      m_ptr_on_device  = rhs.m_ptr_on_device ;
      memory_space::increment( m_ptr_on_device );

      return *this ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool () const
  { return 0 != m_ptr_on_device ; }

  /** \brief  Query if view to same memory */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const View & rhs ) const
  {
    return m_ptr_on_device == rhs.m_ptr_on_device && m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const View & rhs ) const
  {
    return m_ptr_on_device != rhs.m_ptr_on_device || m_shape != rhs.m_shape ;
  }

  /** \brief  Query if view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const View<rhsDataType,rhsLayout,rhsMemory> & rhs ) const
  {
    return m_ptr_on_device == rhs.m_ptr_on_device && m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const View<rhsDataType,rhsLayout,rhsMemory> & rhs ) const
  {
    return m_ptr_on_device != rhs.m_ptr_on_device || m_shape != rhs.m_shape ;
  }

  /*------------------------------------------------------------------*/

private:

  View( value_type * ptr )
    : m_shape()
    , m_ptr_on_device( ptr )
    { memory_space::increment( m_ptr_on_device ); }

  shape_type   m_shape ;
  value_type * m_ptr_on_device ;

  template< class , class , class > friend class View ;
  template< class Dst , class Src >  friend class Impl::Factory ;
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

