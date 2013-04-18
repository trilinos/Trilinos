/*
//@HEADER
// ************************************************************************
//
//                             KokkosArray
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_OPENMP_HPP
#define KOKKOSARRAY_OPENMP_HPP

#include <omp.h>
#include <cstddef>
#include <KokkosArray_Host.hpp>
#include <KokkosArray_Layout.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {

/// \class OpenMP
/// \brief KokkosArray device for multicore processors in the host memory space.
class OpenMP {
public:
  //! \name Type declarations that all KokkosArray devices must provide.
  //@{

  typedef OpenMP                type ;
  typedef OpenMP                device_type ;
  typedef OpenMP                layout_type ;
  typedef HostSpace::size_type  size_type ;
  typedef HostSpace             memory_space ;
  typedef LayoutRight           array_layout ;

  //@}
  //! \name Functions that all KokkosArray devices must implement.
  //@{

  /** \brief  Set the device in a "sleep" state. A noop for OpenMP.  */
  static bool sleep();

  /** \brief Wake the device from the 'sleep' state. A noop for OpenMP. */
  static bool wake();

  /** \brief Wait until all dispatched functors complete. A noop for OpenMP. */
  static void fence() {}

  /// \brief Free any resources being consumed by the device.
  static void finalize();

  /** \brief  Initialize the device.
   *
   *  1) If the hardware locality library is enabled and OpenMP has not
   *     already bound threads then bind OpenMP threads to maximize
   *     core utilization and group for memory hierarchy locality.
   *
   *  2) Allocate a HostThread for each OpenMP thread to hold its
   *     topology and fan in/out data.
   */
  static void initialize( const unsigned gang_count = 1 ,
                          const unsigned worker_per_gang = 0 );

  //------------------------------------

  static void resize_reduce_scratch( unsigned );

  static void * root_reduce_scratch();

  static void assert_ready( const char * const );

  inline static
  Impl::HostThread * get_host_thread()
    { return m_host_threads[ omp_get_thread_num() ]; }

  inline static
  Impl::HostThread * get_host_thread( const unsigned i )
    { return m_host_threads[ i ]; }

private:
  static Impl::HostThread * m_host_threads[ Impl::HostThread::max_thread_count ];
};

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray

#include <OpenMP/KokkosArray_OpenMP_Parallel.hpp>

#endif /* #define KOKKOSARRAY_OPENMP_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

