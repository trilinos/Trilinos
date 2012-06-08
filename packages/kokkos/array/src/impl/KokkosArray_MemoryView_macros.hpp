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

#error "Including <impl/KokkosArray_MemoryView_macros.hpp> without macros defined"

#else

namespace KokkosArray {
namespace Impl {

// Partial specialization of MemoryView, using the memory_space
// typedef provided by the device.
template< typename ValueType >
class MemoryView< ValueType , KOKKOS_MACRO_DEVICE::memory_space > {
public:
  typedef ValueType                          value_type ;
  typedef KOKKOS_MACRO_DEVICE::memory_space  memory_space ;
  typedef KOKKOS_MACRO_DEVICE::size_type     size_type ;
  typedef MemoryView< ValueType , Host >     HostMirror ;

private:

  typedef MemoryManager< memory_space >  memory_manager ;
  typedef memory_manager::view_tracker   view_tracker ;

  ValueType  * m_ptr_on_device ;
  view_tracker m_tracker ;

public:

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

  /** \brief  If the RHS is a different memory space then not equal */
  template< class MemoryRHS >
  inline
  bool operator == ( const MemoryView< ValueType , MemoryRHS > & ) const
  { return false ; }

  /** \brief  If the RHS is a different memory space then not equal */
  template< class MemoryRHS >
  inline
  bool operator != ( const MemoryView< ValueType , MemoryRHS > & ) const
  { return true ; }

  /*------------------------------------------------------------------*/
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~MemoryView()
  {
    memory_manager::clear( m_tracker , m_ptr_on_device );
    m_ptr_on_device = 0 ;
  }

  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView() : m_ptr_on_device(0)
    { memory_manager::init( m_tracker ); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView( const MemoryView & rhs )
    : m_ptr_on_device(0)
    {
      memory_manager::init(   m_tracker );
      memory_manager::assign( m_tracker , rhs.m_tracker );
      m_ptr_on_device = rhs.m_ptr_on_device ;
    }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView & operator = ( const MemoryView & rhs )
    {
      memory_manager::clear(  m_tracker , m_ptr_on_device );
      memory_manager::assign( m_tracker , rhs.m_tracker );
      m_ptr_on_device = rhs.m_ptr_on_device ;
      return *this ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Allocation occurs only on the host */
  inline
  void allocate( size_type count , const std::string & label )
    {
      memory_manager::clear( m_tracker , m_ptr_on_device );

      if ( count ) {
        memory_manager::assign( m_tracker , m_tracker );

        m_ptr_on_device = (value_type *)
          memory_manager::allocate(
            label , typeid(value_type) , sizeof(value_type) , count );
      }
      else {
        m_ptr_on_device = 0 ;
      }
    }

  /*------------------------------------------------------------------*/
  /** \brief  On the host for testing purposes only
   *          can get a count of number of views on the host.
   */
  size_t test_support_view_count() const
    { return m_tracker.test_support_view_count(); }
};

//----------------------------------------------------------------------------

template<>
class MemoryView< void , KOKKOS_MACRO_DEVICE::memory_space > {
public:
  typedef void                               value_type ;
  typedef KOKKOS_MACRO_DEVICE::memory_space  memory_space ;
  typedef KOKKOS_MACRO_DEVICE::size_type     size_type ;
  typedef MemoryView< void , Host >          HostMirror ;

  /*------------------------------------------------------------------*/

#if defined(KOKKOS_MACRO_DEVICE_FUNCTION)

  /** \brief  Query value at offset */
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator[]( const iType & ) const {}

#endif /* defined(KOKKOS_MACRO_DEVICE_FUNCTION) */

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void * ptr_on_device() const { return 0 ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool() const { return false ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const MemoryView & rhs ) const
    { return true ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const MemoryView & rhs ) const
    { return false ; }

  /** \brief  If the RHS is a different memory space then not equal */
  template< class MemoryRHS >
  inline
  bool operator == ( const MemoryView< void , MemoryRHS > & ) const
    { return false ; }

  /** \brief  If the RHS is a different memory space then not equal */
  template< class MemoryRHS >
  inline
  bool operator != ( const MemoryView< void , MemoryRHS > & ) const
    { return true ; }

  /*------------------------------------------------------------------*/
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~MemoryView() {}

  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView() {}

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView( const MemoryView & ) {}

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MemoryView & operator = ( const MemoryView & ) { return *this ; }

  /*------------------------------------------------------------------*/
  /** \brief  Allocation occurs only on the host */
  inline
  void allocate( size_type count , const std::string & label ) {}

  /*------------------------------------------------------------------*/
  /** \brief  On the host for testing purposes only
   *          can get a count of number of views on the host.
   */
  size_t test_support_view_count() const { return 0 ; }
};

} // namespace Impl
} // namespace KokkosArray

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif

