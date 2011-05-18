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

#ifndef KOKKOS_DEVICECUDA_HPP
#define KOKKOS_DEVICECUDA_HPP

#include <iosfwd>
#include <typeinfo>
#include <Kokkos_MemoryView.hpp>
#include <Kokkos_MDArrayIndexMaps.hpp>
#include <impl/Kokkos_ViewTracker.hpp>

#define KOKKOS_DEVICE_CUDA  Kokkos::DeviceCuda

#include <impl/Kokkos_DeviceCuda_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class DeviceCuda {
private:

  static void * allocate_memory( size_t bytes , const std::string & label ,
                                 const std::type_info & type );

  static void deallocate_memory( void * );

  static unsigned m_launching_kernel ;

public:

  /** \brief  On the cuda device use unsigned int for indexing */
  typedef unsigned int         size_type ;
  typedef MDArrayIndexMapLeft  default_mdarray_map ;

  /*--------------------------------*/
  /** \brief  Clear the memory view setting it to the NULL view.
   *          If this is the last view to this allocated memory
   *          then deallocate this allocated memory.
   */
  template< typename ValueType >
  static
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void clear_memory_view( MemoryView< ValueType , DeviceCuda > & lhs )
    {
#if ! defined( KOKKOS_MACRO_DEVICE_FUNCTION )
      // Memory management only available on the host side.
      // If compiling for the device then omit memory management.
      if ( lhs.m_tracker.remove_and_query_is_last() ) {
        deallocate_memory( lhs.m_ptr_on_device );
      }
#endif
      lhs.m_ptr_on_device = 0 ;
    }

  /** \brief  Assign the 'lhs' view to be another view of the 'rhs' view.
   *          Clear the 'lhs' view before the assignment.
   */
  template< typename ValueType >
  static
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void assign_memory_view(       MemoryView< ValueType , DeviceCuda > & lhs ,
                           const MemoryView< ValueType , DeviceCuda > & rhs )
    {
      clear_memory_view( lhs );
#if ! defined( KOKKOS_MACRO_DEVICE_FUNCTION )
      // Memory management only available on the host side.
      // If compiling for the device then omit memory management.
      // If launching a kernel then the view is untracked.
      if ( ! m_launching_kernel ) {
        lhs.m_tracker.insert( rhs.m_tracker );
      }
#endif
      lhs.m_ptr_on_device = rhs.m_ptr_on_device ;
    }

  /** \brief  Allocate memory to be viewed by 'lhs' */
  template< typename ValueType >
  static
  void allocate_memory_view( MemoryView< ValueType , DeviceCuda > & lhs ,
                             size_t count , const std::string & label )
    {
      clear_memory_view( lhs );  
      lhs.m_ptr_on_device = allocate_memory( sizeof(ValueType) * count , label , typeid(ValueType) );
      lhs.m_tracker.insert( lhs );
    }

  /** \brief  Print information about allocate memory */
  void print_memory_view( std::ostream & ) const ;

  /*--------------------------------*/
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <impl/Kokkos_MemoryView_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>

#endif /* #ifndef KOKKOS_DEVICECUDA_HPP */

