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

#error "Including <Kokkos_CrsArray_macros.hpp> without macros defined"

#else

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename SizeType >
class CrsArray< void , KOKKOS_MACRO_DEVICE , SizeType >
{
public:
  typedef KOKKOS_MACRO_DEVICE  device_type ;
  typedef SizeType             size_type ;
  typedef void                 array_type ;
  typedef void                 value_type ;

  typedef CrsArray< void ,
                    HostMapped< device_type >::type ,
                    size_type > HostMirror ;

  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type row_count() const { return m_row_count ; }

  template< typename iTypeRow >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type row_entry_begin( const iTypeRow & row ) const
    { return m_row_map[row]; }

  template< typename iTypeRow >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type row_entry_end( const iTypeRow & row ) const
    { return m_row_map[row+1]; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const size_type * row_map() const
    { return m_row_map.ptr_on_device(); }
  
  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_rank() const { return 0 ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_dimension( const iType & ) const { return 0 ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const void * data() const { return 0 ; }
  
  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray()
  : m_row_map()
  , m_row_count(0)
  {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray( const CrsArray & rhs )
  : m_row_map(     rhs.m_row_map )
  , m_row_count(   rhs.m_row_count )
  {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray & operator = ( const CrsArray & rhs )
  {
    m_row_map     = rhs.m_row_map ;
    m_row_count   = rhs.m_row_count ;
    return *this ;
  }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~CrsArray() { m_row_count = 0 ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const
  { return m_row_map.operator bool(); }

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsArray & rhs ) const
  { return m_row_map.ptr_on_device() == rhs.m_row_map.ptr_on_device() ; }

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsArray & rhs ) const
  { return m_row_map.ptr_on_device() != rhs.m_row_map.ptr_on_device() ; }

private:

  typedef device_type::memory_space  memory_space ;

  Impl::MemoryView< size_type ,  memory_space > m_row_map ;
  size_type                                     m_row_count ;

  template< typename , class , typename > friend class CrsArray ;
  template< class Dst , class Src >  friend class Impl::Factory ;
};

//----------------------------------------------------------------------------

template< class ArrayType , typename SizeType >
class CrsArray< ArrayType , KOKKOS_MACRO_DEVICE , SizeType >
{
public:
  typedef KOKKOS_MACRO_DEVICE                                  device_type ;
  typedef SizeType                                             size_type ;
  typedef ArrayType                                            array_type ;
  typedef typename Impl::remove_all_extents<array_type>::type  value_type ;

  typedef CrsArray< array_type ,
                    HostMapped< device_type >::type ,
                    size_type > HostMirror ;

  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type row_count() const { return m_row_count ; }

  template< typename iTypeRow >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type row_entry_begin( const iTypeRow & row ) const
    { return m_row_map[row]; }

  template< typename iTypeRow >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type row_entry_end( const iTypeRow & row ) const
    { return m_row_map[row+1]; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const size_type * row_map() const
    { return m_row_map.ptr_on_device(); }
  
  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_rank() const { return m_index_map.rank(); }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_dimension( const iType & rank ) const
    { return m_index_map.dimension( rank ); }

  template< typename iTypeEntry >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeEntry & entry ) const
    { return m_data[ entry ]; }

  template< typename iTypeEntry , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeEntry & entry ,
                           const iType1     & i1 ) const
    { return m_data[ m_index_map.offset( entry , i1 ) ]; }

  template< typename iTypeEntry , typename iType1 , typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeEntry & entry ,
                           const iType1     & i1 ,
                           const iType2     & i2 ) const
    { return m_data[ m_index_map.offset( entry , i1 , i2 ) ]; }

  template< typename iTypeEntry , typename iType1 , typename iType2 ,
            typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeEntry & entry ,
                           const iType1     & i1 ,
                           const iType2     & i2 ,
                           const iType3     & i3 ) const
    { return m_data[ m_index_map.offset( entry , i1 , i2 , i3 ) ]; }

  template< typename iTypeEntry , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeEntry & entry ,
                           const iType1     & i1 ,
                           const iType2     & i2 ,
                           const iType3     & i3 ,
                           const iType4     & i4 ) const
    { return m_data[ m_index_map.offset( entry , i1 , i2 , i3 , i4 ) ]; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const value_type * data() const
    { return m_data.ptr_on_device(); }
  
  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray()
  : m_row_map()
  , m_data()
  , m_index_map()
  , m_row_count(0)
  {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray( const CrsArray & rhs )
  : m_row_map(     rhs.m_row_map )
  , m_data(        rhs.m_data )
  , m_index_map(   rhs.m_index_map )
  , m_row_count(   rhs.m_row_count )
  {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray & operator = ( const CrsArray & rhs )
  {
    m_row_map     = rhs.m_row_map ;
    m_data        = rhs.m_data ;
    m_index_map   = rhs.m_index_map ;
    m_row_count   = rhs.m_row_count ;
    return *this ;
  }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~CrsArray() { m_row_count = 0 ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const
  { return m_row_map.operator bool(); }

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsArray & rhs ) const
  { return m_row_map.ptr_on_device() == rhs.m_row_map.ptr_on_device() ; }

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsArray & rhs ) const
  { return m_row_map.ptr_on_device() != rhs.m_row_map.ptr_on_device() ; }

  /*------------------------------------------------------------------*/

  template< class DevRHS >
  bool operator == ( const CrsArray< ArrayType , DevRHS , SizeType > & rhs ) const
  {
    return Impl::SameType< device_type , DevRHS >::value &&
           m_row_map.ptr_on_device() == rhs.m_row_map.ptr_on_device() ;
  }

  template< class DevRHS >
  bool operator != ( const CrsArray< ArrayType , DevRHS , SizeType > & rhs ) const
  { return ! operator == ( rhs ); }

  /*------------------------------------------------------------------*/

private:

  typedef device_type::memory_space  memory_space ;

  static const unsigned Rank = Impl::rank< array_type >::value ;

  typedef typename device_type::IndexMap<
    Rank + 1 ,
    0 < Rank ? Impl::extent< array_type, 0 >::value : 0 ,
    1 < Rank ? Impl::extent< array_type, 1 >::value : 0 ,
    2 < Rank ? Impl::extent< array_type, 2 >::value : 0 ,
    3 < Rank ? Impl::extent< array_type, 3 >::value : 0 ,
    4 < Rank ? Impl::extent< array_type, 4 >::value : 0 ,
    5 < Rank ? Impl::extent< array_type, 5 >::value : 0 ,
    6 < Rank ? Impl::extent< array_type, 6 >::value : 0 >::type
  internal_index_map_type ;

  Impl::MemoryView< size_type ,  memory_space > m_row_map ;
  Impl::MemoryView< value_type,  memory_space > m_data ;
  internal_index_map_type                       m_index_map ;
  size_type                                     m_row_count ;

  template< typename , class , typename > friend class CrsArray ;
  template< class Dst , class Src >  friend class Impl::Factory ;
};


//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

