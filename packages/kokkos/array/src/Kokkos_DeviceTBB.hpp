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

#ifndef KOKKOS_DEVICETBB_HPP
#define KOKKOS_DEVICETBB_HPP

#include <Kokkos_DeviceHost.hpp>

#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>

#define KOKKOS_DEVICE_TBB  Kokkos::DeviceTBB

/*--------------------------------------------------------------------------*/

namespace Kokkos {
	
class MDArrayIndexMapRight;

class DeviceTBB {
private:
 
  static unsigned m_launching_kernel ;

public:

  /** \brief  On the TBB device use size_t for indexing */
  typedef size_t               size_type ;

  /** \brief  The TBB device uses the Host memory space */
  typedef DeviceHost           memory_space ;

  /** \brief  Default mdarray map is index from right */
  typedef MDArrayIndexMapRight  default_mdarray_map ;

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

#endif /* #define KOKKOS_DEVICETBB_HPP */

