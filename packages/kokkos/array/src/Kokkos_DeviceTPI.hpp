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

#ifndef KOKKOS_DEVICETPI_HPP
#define KOKKOS_DEVICETPI_HPP

#include <Kokkos_DeviceHost.hpp>

#define KOKKOS_DEVICE_TPI  Kokkos::DeviceTPI

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class DeviceTPI {
private:

  static unsigned m_launching_kernel ;

public:

  /** \brief  On the TPI device use size_t for indexing */
  typedef size_t                size_type ;

  /** \brief  The TPI device uses the Host memory space */
  typedef HostMemory            memory_space ;

  /** \brief  Default mdarray map is index from right */
  typedef Impl::MDArrayIndexMapRight  mdarray_map ;

  /*--------------------------------*/

  static void initialize( size_type nthreads );

  static void finalize();

  /*--------------------------------*/

  static void set_dispatch_functor();
  static void clear_dispatch_functor();
  static void wait_functor_completion() {}

  /*--------------------------------*/
};

} // namespace Kokkos

#include <Kokkos_DeviceTPI_macros.hpp>
#include <impl/Kokkos_BasicFunctors_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

#endif /* #define KOKKOS_DEVICETPI_HPP */

