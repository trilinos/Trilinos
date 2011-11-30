/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Kokkos_CrsArray_macros.hpp> without macros defined"

#else

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Compressed row storage array.
 *
 *  Map ( I , J ) -> K
 *  where I in [ 0 .. N ) and J in [ 0 .. M_I )
 */
template< typename ValueType >
class CrsArray< ValueType , KOKKOS_MACRO_DEVICE > {
public:
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef ValueType               value_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsArray< value_type , Host > HostView ;

  /*------------------------------------------------------------------*/
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type row_dimension() const { return m_row_count ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const size_type * offset_array() const
    { return m_offset.ptr_on_device(); }
  
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * values_array() const
    { return m_values.ptr_on_device(); }
  
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type size() const
    { return offset_array()[ m_row_count ]; }

  /*------------------------------------------------------------------*/
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type column_dimension( const iType & row ) const
    {
      const size_type * const offset = offset_array() + row ;
      return offset[1] - offset[0] ;
    }

  template< typename iType , typename jType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iType & row , const jType & j ) const
    {
      const size_type offset = offset_array()[row] + j ;
      return values_array()[ offset ];
    }

  typedef std::pair<value_type*,value_type*> value_range_type ;

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_range_type operator()( const iType & row ) const
    {
      const size_type  * const offset = offset_array() + row ;
      const value_type * const values = values_array();
      return value_range_type( values + offset[0] , values + offset[1] );
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray()
  : m_offset()
  , m_values()
  , m_row_count(0)
  {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray( const CrsArray & rhs )
  : m_offset(    rhs.m_offset )
  , m_values(    rhs.m_values )
  , m_row_count( rhs.m_row_count )
  {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsArray & operator = ( const CrsArray & rhs )
  {
    m_offset    = rhs.m_offset ;
    m_values    = rhs.m_values ;
    m_row_count = rhs.m_row_count ;
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
  { return m_offset.operator bool(); }

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsArray & rhs ) const
  { return m_offset.operator==( rhs.m_offset ); }

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsArray & rhs ) const
  { return m_offset.operator!=( rhs.m_offset ); }

private:

  typedef device_type::memory_space  memory_space ;

  Impl::MemoryView< size_type ,  memory_space > m_offset ;
  Impl::MemoryView< value_type , memory_space > m_values ;
  size_type                                     m_row_count ;

  template< class D , typename V > friend class Impl::CreateCrsArray ;
  template< class Dst , class Src > friend class Impl::DeepCopy ;
};

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

