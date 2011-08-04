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

#ifndef KOKKOS_DEVICEFERRY_HPP
#define KOKKOS_DEVICEFERRY_HPP

#include <iosfwd>
#include <typeinfo>
#include <vector>

#include <Kokkos_MemoryView.hpp>
#include <impl/Kokkos_ViewTracker.hpp>

#define KOKKOS_DEVICE_FERRY  Kokkos::DeviceFerry

#include <Kokkos_DeviceFerry_macros.hpp>


/*--------------------------------------------------------------------------*/

namespace Kokkos {

class DeviceFerry {
private:

  static void load_parallel_functor( const void * , size_t );

  static const void * get_parallel_functor();

	
  static void * allocate_memory( const std::string & label ,
                                 const std::type_info & type ,
                                 const size_t member_size ,
                                 const size_t member_count );

  static void deallocate_memory( void * );

  static unsigned m_launching_kernel ;

public:

  /** \brief  On the ferry device use unsigned int for indexing */
  typedef unsigned int         size_type ;

  /** \brief  The Ferry device uses the Ferry memory space */
  typedef DeviceFerry          memory_space ;

  /** \brief  Default mdarray map is index from right */
  typedef Impl::MDArrayIndexMapRight  mdarray_map ;

  /*--------------------------------*/

  static void initialize( );
  
  /*--------------------------------*/
  /** \brief  Clear the memory view setting it to the NULL view.
   *          If this is the last view to this allocated memory
   *          then deallocate this allocated memory.
   */
  template< typename ValueType >
  static void clear_memory_view( MemoryView< ValueType , DeviceFerry > & lhs )
    {
#if ! defined( __MIC__ )
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
  static void assign_memory_view(       MemoryView< ValueType , DeviceFerry > & lhs ,
                           const MemoryView< ValueType , DeviceFerry > & rhs )
    {
      clear_memory_view( lhs );
#if ! defined( __MIC__ )
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
  void allocate_memory_view( MemoryView< ValueType , DeviceFerry > & lhs ,
                             size_t count , const std::string & label )
    {
      clear_memory_view( lhs );  
      lhs.m_ptr_on_device = (ValueType *)
        allocate_memory( label, typeid(ValueType), sizeof(ValueType), count );
      lhs.m_tracker.insert( lhs.m_tracker );
    }

  /** \brief  Print information about allocate memory */
  static void print_memory_view( std::ostream & );

  /*--------------------------------*/

  static void set_dispatch_functor();
  static void clear_dispatch_functor();
  static void wait_functor_completion();

  /*--------------------------------*/
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <impl/Kokkos_MemoryView_macros.hpp>

#include <Kokkos_DeviceClear_macros.hpp>


#endif /* #ifndef KOKKOS_DEVICEFERRY_HPP */

