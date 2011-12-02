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

#error "Including <Kokkos_CrsMap_macros.hpp> without macros defined"

#else

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename SizeType >
class CrsMap< KOKKOS_MACRO_DEVICE , SizeType > {
public:
  typedef KOKKOS_MACRO_DEVICE  device_type ;
  typedef SizeType             size_type ;

  /*------------------------------------------------------------------*/
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type first_count() const { return m_first_count ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type size() const { return m_total_count ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const size_type * offset_array() const
    { return m_memory.ptr_on_device(); }
  
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type second_count( const iType & first ) const
    {
      const size_type * const offset = m_memory.ptr_on_device() + first ;
      return offset[1] - offset[0] ;
    }

  template< typename iType , typename jType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type operator()( const iType & first , const jType & second ) const
    { return m_memory.ptr_on_device()[first] + second ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  std::pair<size_type,size_type> range( const iType & first ) const
    {
      const size_type * const offset = m_memory.ptr_on_device() + first ;
      return std::pair<size_type,size_type>( offset[0] , offset[1] );
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsMap()
  : m_memory()
  , m_first_count(0)
  , m_total_count(0)
  {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsMap( const CrsMap & rhs )
  : m_memory(    rhs.m_memory )
  , m_first_count( rhs.m_first_count )
  , m_total_count( rhs.m_total_count )
  {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  CrsMap & operator = ( const CrsMap & rhs )
  {
    m_memory    = rhs.m_memory ;
    m_first_count = rhs.m_first_count ;
    m_total_count = rhs.m_total_count ;
    return *this ;
  }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~CrsMap() { m_first_count = 0 ; m_total_count = 0 ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const
  { return m_memory.operator bool(); }

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsMap & rhs ) const
  { return m_memory.ptr_on_device() == rhs.m_memory.ptr_on_device() ; }

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsMap & rhs ) const
  { return m_memory.ptr_on_device() != rhs.m_memory.ptr_on_device() ; }

  /** \brief  Query if view to same memory */
  template< class DeviceRHS , typename SizeTypeRHS >
  bool operator == ( const CrsMap< DeviceRHS , SizeTypeRHS > & rhs ) const
    {
      return Impl::SameType< typename device_type::memory_space ,
                             typename DeviceRHS  ::memory_space >::value &&
             Impl::SameType< size_type , SizeTypeRHS >::value &&
             m_memory.ptr_on_device() ==
               (size_type *) rhs.m_memory.ptr_on_device();
    }

  /** \brief  Query if not view to same memory */
  template< class DeviceRHS , typename SizeTypeRHS >
  bool operator !=
    ( const CrsMap< DeviceRHS , SizeTypeRHS > & rhs ) const
    { return ! operator==( rhs ); }

private:

  typedef device_type::memory_space  memory_space ;

  Impl::MemoryView< size_type ,  memory_space > m_memory ;
  size_type                                     m_first_count ;
  size_type                                     m_total_count ;

  template< class Dev , typename I > friend class CrsMap ;
  template< class Dst , class Src >  friend class Impl::CreateCrsMap ;
};

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

