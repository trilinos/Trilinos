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

#error "Including <impl/Kokkos_MemoryView_macros.hpp> without macros defined"

#else

namespace Kokkos {

template< typename ValueType >
class MemoryView< ValueType , KOKKOS_MACRO_DEVICE > {
private:

  ValueType       * m_ptr_on_device ;
  Impl::ViewTracker m_tracker ;

  friend class KOKKOS_MACRO_DEVICE ;

  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView( const MemoryView & rhs );

  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView & operator = ( const MemoryView & rhs );

public:

  typedef ValueType           value_type ;
  typedef KOKKOS_MACRO_DEVICE device_type ;

  /*------------------------------------------------------------------*/

#if defined(KOKKOS_MACRO_DEVICE_FUNCTION)

  /** \brief  Query value at offset */
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator[]( const iType & i ) const
  { return m_ptr_on_device[ i ]; }

#endif /* defined(KOKKOS_MACRO_DEVICE_FUNCTION) */

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * ptr_on_device() const
  { return m_ptr_on_device ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool() const
  { return 0 != m_ptr_on_device ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const MemoryView & rhs ) const
  { return m_ptr_on_device == rhs.m_ptr_on_device ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const MemoryView & rhs ) const
  { return m_ptr_on_device != rhs.m_ptr_on_device ; }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView() : m_ptr_on_device(0) { m_tracker.next = 0 ; }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~MemoryView()
  { device_type::clear_memory_view( *this ); }

  /*------------------------------------------------------------------*/
  /** \brief  On the host for testing purposes only
   *          can get a count of number of views on the host.
   */
  size_t test_support_view_count() const
    { return m_tracker.test_support_view_count(); }
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#endif

