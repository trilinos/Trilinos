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

#error "Including <KokkosArray_Array_macros.hpp> without macros defined"

#else

namespace KokkosArray {

template< class ArrayType >
class Array< ArrayType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE                                  device_type ;
  typedef device_type::size_type                               size_type ;
  typedef ArrayType                                            array_type ;
  typedef typename Impl::remove_all_extents<array_type>::type  value_type ;

  typedef Array< array_type ,
                    HostMapped< device_type >::type > HostMirror ;

  /*------------------------------------------------------------------*/

  static const unsigned Rank = Impl::rank<array_type>::value + 1 ;

  template< unsigned I >
  struct Dimension {
    static const unsigned value = Impl::extent<array_type,I-1>::value ;
  };

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return m_index_map.rank(); }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & rank ) const
    { return m_index_map.dimension( rank ); }

  /** \brief  Because memory is contiguous this is exposed */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * ptr_on_device() const
  {
    // TBD: If memory is not contiguous and can throw then throw !
    return m_data.ptr_on_device();
  }

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  template< typename iTypeP >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ) const
    { return m_data[ iP ]; }

  template< typename iTypeP , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ) const
    { return m_data[ m_index_map.offset( iP , i1 ) ]; }

  template< typename iTypeP , typename iType1 , typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ,
                           const iType2 & i2 ) const
    { return m_data[ m_index_map.offset( iP , i1 , i2 ) ]; }

  template< typename iTypeP , typename iType1 , typename iType2 ,
            typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ) const
    { return m_data[ m_index_map.offset( iP , i1 , i2 , i3 ) ]; }

  template< typename iTypeP , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ) const
    { return m_data[ m_index_map.offset( iP , i1 , i2 , i3 , i4 ) ]; }

  template< typename iTypeP , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ) const
    { return m_data[ m_index_map.offset( iP , i1 , i2 , i3 , i4 , i5 ) ]; }

  template< typename iTypeP , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ) const
    { return m_data[ m_index_map.offset( iP , i1 , i2 , i3 , i4 , i5 , i6 ) ]; }

  template< typename iTypeP , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ,
                           const iType6 & i7 ) const
    { return m_data[ m_index_map.offset( iP, i1, i2, i3, i4, i5, i6, i7 ) ]; }

#endif /* #if defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  Array()
    : m_data()
    , m_index_map()
    {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  Array( const Array & rhs )
    : m_data(        rhs.m_data )
    , m_index_map(   rhs.m_index_map )
    {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  Array & operator = ( const Array & rhs )
  {
    m_data        = rhs.m_data ;
    m_index_map   = rhs.m_index_map ;
    return *this ;
  }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~Array() {}

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const
  { return m_data.operator bool(); }

  /** \brief  Query if view to same memory */
  bool operator == ( const Array & rhs ) const
  { return m_data.ptr_on_device() == rhs.m_data.ptr_on_device() ; }

  /** \brief  Query if not view to same memory */
  bool operator != ( const Array & rhs ) const
  { return m_data.ptr_on_device() != rhs.m_data.ptr_on_device() ; }

  /*------------------------------------------------------------------*/

  template< class DevRHS >
  bool operator == ( const Array< ArrayType , DevRHS > & rhs ) const
  {
    return Impl::SameType< device_type , DevRHS >::value &&
           m_data.ptr_on_device() == rhs.m_data.ptr_on_device() ;
  }

  template< class DevRHS >
  bool operator != ( const Array< ArrayType , DevRHS > & rhs ) const
  { return ! operator == ( rhs ); }

  /*------------------------------------------------------------------*/

private:

  typedef device_type::memory_space  memory_space ;

  typedef typename device_type::IndexMap<
    Rank ,
    1 < Rank ? Impl::extent< array_type, 0 >::value : 0 ,
    2 < Rank ? Impl::extent< array_type, 1 >::value : 0 ,
    3 < Rank ? Impl::extent< array_type, 2 >::value : 0 ,
    4 < Rank ? Impl::extent< array_type, 3 >::value : 0 ,
    5 < Rank ? Impl::extent< array_type, 4 >::value : 0 ,
    6 < Rank ? Impl::extent< array_type, 5 >::value : 0 ,
    7 < Rank ? Impl::extent< array_type, 6 >::value : 0 >::type
  internal_index_map_type ;

  Impl::MemoryView< value_type,  memory_space > m_data ;
  internal_index_map_type                       m_index_map ;

  template< typename , class > friend class Array ;
  template< class Dst , class Src >  friend class Impl::Factory ;
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

