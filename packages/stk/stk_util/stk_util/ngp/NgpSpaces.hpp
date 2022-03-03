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

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif

namespace stk {
namespace ngp {

#ifdef KOKKOS_ENABLE_CUDA
using ExecSpace = Kokkos::Cuda;
#elif defined(KOKKOS_ENABLE_HIP)
using ExecSpace = Kokkos::Experimental::HIP;
#elif defined(KOKKOS_ENABLE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
#else
using ExecSpace = Kokkos::Serial;
#endif

#ifdef KOKKOS_ENABLE_CUDA
using HostExecSpace = Kokkos::Serial;
#elif defined(KOKKOS_ENABLE_HIP)
using HostExecSpace = Kokkos::Serial;
#elif defined(KOKKOS_ENABLE_OPENMP)
using HostExecSpace = Kokkos::OpenMP;
#else
using HostExecSpace = Kokkos::Serial;
#endif

#ifdef KOKKOS_ENABLE_CUDA
using MemSpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ENABLE_HIP)
using MemSpace = Kokkos::Experimental::HIPSpace;
#elif defined(KOKKOS_ENABLE_OPENMP)
using MemSpace = Kokkos::OpenMP;
#else
using MemSpace = Kokkos::HostSpace;
#endif

#ifdef KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_CUDA_UVM
using UVMMemSpace = Kokkos::CudaUVMSpace;
#else
using UVMMemSpace = Kokkos::CudaHostPinnedSpace;
#endif
#elif defined(KOKKOS_ENABLE_HIP)
using UVMMemSpace = Kokkos::Experimental::HIPHostPinnedSpace;
#elif defined(KOKKOS_ENABLE_OPENMP)
using UVMMemSpace = Kokkos::OpenMP;
#else
using UVMMemSpace = Kokkos::HostSpace;
#endif

#ifdef KOKKOS_ENABLE_CUDA
using HostPinnedSpace = Kokkos::CudaHostPinnedSpace;
#elif defined(KOKKOS_ENABLE_HIP)
using HostPinnedSpace = Kokkos::Experimental::HIPHostPinnedSpace;
#else
using HostPinnedSpace = MemSpace;
#endif

using ScheduleType = Kokkos::Schedule<Kokkos::Dynamic>;

}
}

#endif
