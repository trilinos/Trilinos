// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_NGP_NGPSPACES_HPP
#define STK_NGP_NGPSPACES_HPP

#include "stk_util/stk_config.h"
#include <Kokkos_Core.hpp>

namespace stk {
namespace ngp {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

#ifndef KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
#ifndef _MSC_VER
#warning "Kokkos::SharedHostPinnedSpace is not defined."
#else
#pragma message("Kokkos::SharedHostPinnedSpace is not defined.")
#endif
using HostPinnedSpace = Kokkos::HostSpace;
#else
using HostPinnedSpace = Kokkos::SharedHostPinnedSpace;
#endif

#ifndef KOKKOS_HAS_SHARED_SPACE
#ifndef _MSC_VER
#warning "Kokkos::SharedSpace is not defined."
#else
#pragma message("Kokkos::SharedSpace is not defined.")
#endif
using UVMMemSpace = Kokkos::HostSpace;
#else
using UVMMemSpace = Kokkos::SharedSpace;
#endif

#ifdef KOKKOS_ENABLE_CUDA
using MemSpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ENABLE_HIP)
using MemSpace = Kokkos::HIPSpace;
#else
using MemSpace = ExecSpace::memory_space;
#endif

using HostMemSpace = HostExecSpace::memory_space;

#ifdef KOKKOS_ENABLE_HIP
template <typename ExecutionSpace>
using RangePolicy = Kokkos::RangePolicy<ExecutionSpace, Kokkos::LaunchBounds<128, 1>>;

using HostRangePolicy = Kokkos::RangePolicy<HostExecSpace, Kokkos::LaunchBounds<128, 1>>;
using DeviceRangePolicy = Kokkos::RangePolicy<ExecSpace, Kokkos::LaunchBounds<128, 1>>;

template <typename ExecutionSpace>
using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace, Kokkos::LaunchBounds<128, 1>>;

using HostTeamPolicy = Kokkos::TeamPolicy<HostExecSpace, Kokkos::LaunchBounds<128, 1>>;
using DeviceTeamPolicy = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<128, 1>>;
#else
template <typename ExecutionSpace>
using RangePolicy = Kokkos::RangePolicy<ExecutionSpace>;

using HostRangePolicy = Kokkos::RangePolicy<HostExecSpace>;
using DeviceRangePolicy = Kokkos::RangePolicy<ExecSpace>;

template <typename ExecutionSpace>
using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;

using HostTeamPolicy = Kokkos::TeamPolicy<HostExecSpace>;
using DeviceTeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
#endif

using ScheduleType = Kokkos::Schedule<Kokkos::Dynamic>;


// Detect if the host and device memory spaces are unified
constexpr bool DeviceAccessibleFromHost = Kokkos::SpaceAccessibility<HostExecSpace, ExecSpace::memory_space>::accessible;


struct HostSpace {
  using exec_space = HostExecSpace;
  using mem_space = HostMemSpace;
};


#if defined(STK_ENABLE_GPU)

  struct DeviceSpace {
    using exec_space = ExecSpace;
    using mem_space = MemSpace;
  };

  struct UVMDeviceSpace {
    using exec_space = ExecSpace;
    using mem_space = UVMMemSpace;
  };

  struct HostPinnedDeviceSpace {
    using exec_space = ExecSpace;
    using mem_space = HostPinnedSpace;
  };

#else

  #if defined(STK_USE_DEVICE_MESH)

    // Fake device space that lives on host, for testing and debugging purposes only (wasteful and poorly-performing)
    struct DeviceSpace {
      using exec_space = HostExecSpace;
      using mem_space = HostMemSpace;
    };

    struct UVMDeviceSpace {
      using exec_space = HostExecSpace;
      using mem_space = HostMemSpace;
    };

    struct HostPinnedDeviceSpace {
      using exec_space = HostExecSpace;
      using mem_space = HostMemSpace;
    };

  #else

    // Alias device to host, to disable specialized device code
    using DeviceSpace = HostSpace;
    using UVMDeviceSpace = HostSpace;
    using HostPinnedDeviceSpace = HostSpace;

  #endif

#endif





}
}

#endif
