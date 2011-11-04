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

#ifndef KOKKOS_DEVICENUMA_HPP
#define KOKKOS_DEVICENUMA_HPP

#include <Kokkos_DeviceHost.hpp>

#define KOKKOS_DEVICE_NUMA  Kokkos::DeviceNUMA

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class DeviceNUMAWorker ;
class DeviceNUMAThread ;
}
}

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class DeviceNUMA {
public:

  /** \brief  On the NUMA device use size_t for indexing */
  typedef size_t                size_type ;

  /** \brief  The NUMA device uses the Host memory space */
  typedef HostMemory            memory_space ;

  /** \brief  Default mdarray map is index from right */
  typedef Impl::MDArrayIndexMapRight  mdarray_map ;

  /*--------------------------------*/

  static void wait_functor_completion() {}

  static void execute( const Impl::DeviceNUMAWorker & );

  /*--------------------------------*/
  enum UtilizationOptions { FULL ///< Use every core of every node
                          , MOST ///< Use all but one core of every node
                          };

  static void initialize( UtilizationOptions );
  static void finalize();

  static void block();
  static void unblock();
};

} // namespace Kokkos

#include <DeviceNUMA/Kokkos_DeviceNUMA_Impl.hpp>
#include <DeviceNUMA/Kokkos_DeviceNUMA_For.hpp>
#include <DeviceNUMA/Kokkos_DeviceNUMA_Reduce.hpp>
#include <DeviceNUMA/Kokkos_DeviceNUMA_views.hpp>

#endif /* #define KOKKOS_DEVICENUMA_HPP */

