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

#if ! defined(KOKKOS_MACRO_IMPL_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_HOST_FUNCTION)           || \
    ! defined(KOKKOS_MACRO_DEVICE_FUNCTION)

#include <macros/Kokkos_Preprocessing_macros.hpp>

#error "Including " ## KOKKOS_MACRO_TO_STRING( __FILE__ ) ## " without macros defined"

#else

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class MemoryView< ValueType , KOKKOS_MACRO_DEVICE > {
private:

  MemoryViewTracker m_tracker ;
  ValueType       * m_ptr_on_device ;

  friend class KOKKOS_MACRO_DEVICE ;

  KOKKOS_MACRO_HOST_FUNCTION
  KOKKOS_MACRO_DEVICE_FUNCTION
  MemoryView( const MemoryView & rhs );

  KOKKOS_MACRO_HOST_FUNCTION
  KOKKOS_MACRO_DEVICE_FUNCTION
  MemoryView & operator = ( const MemoryView & rhs );

public:

  typedef KOKKOS_MACRO_DEVICE device_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value at offset */
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator[]( const iType & i ) const
  { return m_ptr_on_device[ i ]; }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type * ptr_on_device() const
  { return m_ptr_on_device ; }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_HOST_FUNCTION
  KOKKOS_MACRO_DEVICE_FUNCTION
  MemoryView() : m_tracker(), m_ptr_on_device(0) {}

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_HOST_FUNCTION
  KOKKOS_MACRO_DEVICE_FUNCTION
  ~MemoryView()
  { device_type::clear_memory_view( *this ); }
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos

#endif

