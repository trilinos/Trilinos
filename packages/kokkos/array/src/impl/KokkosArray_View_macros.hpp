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
class View< DataType , LayoutSpec , KOKKOS_MACRO_DEVICE::memory_space_new >
{
public:
  typedef DataType             data_type ;
  typedef LayoutSpec           layout_spec ;
  typedef KOKKOS_MACRO_DEVICE  device_type ;

  typedef View< data_type , layout_spec , Host >  HostMirror ;

  typedef typename Impl::remove_all_extents<data_type>::type  value_type ;
  typedef typename LayoutSpec::array_layout                   array_layout ;
  typedef typename device_type::memory_space_new              memory_space ;
  typedef typename device_type::size_type                     size_type ;

private:

  typedef Impl::Shape< typename Impl::remove_const<data_type> ,
                       array_layout > shape_type ;

  typedef Impl::unsigned_< shape_type::rank > Rank ;

  typedef Impl::LayoutMap< shape_type , memory_space > layout_map ;

public:

  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return shape_type::rank ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & rank ) const
    { return Impl::dimension( m_shape , rank ); }

  /** \brief  Because memory is contiguous this is exposed */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * ptr_on_device() const
  { return m_ptr_on_device ; }

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator * () const
    {
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<0> , Rank >::type ok_rank ;
      return *m_ptr_on_device ;
    }

  template< typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ) const
    {
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<1> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0 ) );

      return m_ptr_on_device[ i0 ]; 
    }

  template< typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ) const
    {
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<2> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0, i1 ) );

      return m_ptr_on_device[ layout_map::offset( m_shape, i0, i1 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ) const
    {
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<3> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK( Impl::array_bounds_check( m_shape, i0, i1, i2 ) );

      return m_ptr_on_device[ layout_map::offset( m_shape, i0 , i1 , i2 ) ];
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
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<4> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3 ) );

      return m_ptr_on_device[ layout_map::offset( m_shape, i0, i1, i2, i3 ) ];
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
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<5> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4 ) );

      return m_ptr_on_device[
               layout_map::offset( m_shape, i0, i1, i2, i3, i4 ) ];
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
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<6> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5 ) );

      return m_ptr_on_device[
               layout_map::offset( m_shape, i0, i1, i2, i3, i4, i5 ) ];
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
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<7> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5, i6 ) );

      return m_ptr_on_device[
               layout_map::offset( m_shape, i0, i1, i2, i3, i4, i5, i6 ) ];
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
      typedef typename Impl::StaticAssertSame< Impl::unsigned_<8> , Rank >::type ok_rank ;

      KOKKOS_MACRO_CHECK(
        Impl::array_bounds_check( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 ) );

      return m_ptr_on_device[ layout_map::offset(m_shape,i0,i1,i2,i3,i4,i5,i6,i7) ];
    }

#endif /* #if defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

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
    : m_shape(         rhs.shape )
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

  template< class rhsType , class rhsMapSpec , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View( const View< rhsType , rhsMapSpec , rhsMemory > & rhs )
    : m_shape(         rhs.m_shape ) // Must be same type
    , m_ptr_on_device( rhs.m_ptr_on_device ) // preserves 'const' requirement
    {
      memory_space::increment( m_ptr_on_device );
    }

  template< class rhsType , class rhsMapSpec , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View & operator = ( const View< rhsType , rhsMapSpec , rhsMemory > & rhs )
    {
      memory_space::decrement( m_ptr_on_device );
      m_shape          = rhs.m_shape ; // Must be same type
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
  /** \brief  Query if NULL view */
  operator bool () const
  { return 0 != m_ptr_on_device ; }

  /** \brief  Query if view to same memory */
  bool operator == ( const View & rhs ) const
  { return m_shape == rhs.m_shape && m_ptr_on_device == rhs.ptr_on_device ; }

  /** \brief  Query if not view to same memory */
  bool operator != ( const View & rhs ) const
  { return m_shape != rhs.m_shape || m_ptr_on_device != rhs.ptr_on_device ; }

  /*------------------------------------------------------------------*/

private:

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

