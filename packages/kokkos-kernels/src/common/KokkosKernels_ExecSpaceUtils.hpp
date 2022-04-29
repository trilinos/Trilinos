/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include "Kokkos_Core.hpp"

#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GPU)
#include <level_zero/zes_api.h>
#include <CL/sycl/backend/level_zero.hpp>
#endif

#ifndef _KOKKOSKERNELSUTILSEXECSPACEUTILS_HPP
#define _KOKKOSKERNELSUTILSEXECSPACEUTILS_HPP

namespace KokkosKernels {

namespace Impl {

enum ExecSpaceType {
  Exec_SERIAL,
  Exec_OMP,
  Exec_THREADS,
  Exec_CUDA,
  Exec_HIP,
  Exec_SYCL
};
template <typename ExecutionSpace>
KOKKOS_FORCEINLINE_FUNCTION ExecSpaceType kk_get_exec_space_type() {
  ExecSpaceType exec_space = Exec_SERIAL;
#if defined(KOKKOS_ENABLE_SERIAL)
  if (std::is_same<Kokkos::Serial, ExecutionSpace>::value) {
    exec_space = Exec_SERIAL;
  }
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  if (std::is_same<Kokkos::Threads, ExecutionSpace>::value) {
    exec_space = Exec_THREADS;
  }
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  if (std::is_same<Kokkos::OpenMP, ExecutionSpace>::value) {
    exec_space = Exec_OMP;
  }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<Kokkos::Cuda, ExecutionSpace>::value) {
    exec_space = Exec_CUDA;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  if (std::is_same<Kokkos::Experimental::HIP, ExecutionSpace>::value) {
    exec_space = Exec_HIP;
  }
#endif

#if defined(KOKKOS_ENABLE_SYCL)
  if (std::is_same<Kokkos::Experimental::SYCL, ExecutionSpace>::value) {
    exec_space = Exec_SYCL;
  }
#endif

  return exec_space;
}

////////////////////////////////////////////////////////////////////////////////
// GPU Exec Space Utils
////////////////////////////////////////////////////////////////////////////////

template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION bool kk_is_gpu_exec_space() {
  return false;
}

#ifdef KOKKOS_ENABLE_CUDA
template <>
constexpr KOKKOS_INLINE_FUNCTION bool kk_is_gpu_exec_space<Kokkos::Cuda>() {
  return true;
}
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
kk_is_gpu_exec_space<Kokkos::Experimental::HIP>() {
  return true;
}
#endif

#ifdef KOKKOS_ENABLE_SYCL
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
kk_is_gpu_exec_space<Kokkos::Experimental::SYCL>() {
  return true;
}
#endif

////////////////////////////////////////////////////////////////////////////////
// x86_64 Memory Space Utils
////////////////////////////////////////////////////////////////////////////////

template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION bool kk_is_x86_64_mem_space() {
  return false;
}

#if __x86_64__
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
kk_is_x86_64_mem_space<Kokkos::HostSpace>() {
  return true;
}
#endif  // x86_64 architectures

////////////////////////////////////////////////////////////////////////////////
// A64FX Memory Space Utils
////////////////////////////////////////////////////////////////////////////////

template <typename ExecutionSpace>
constexpr KOKKOS_INLINE_FUNCTION bool kk_is_a64fx_mem_space() {
  return false;
}

#if defined(__ARM_ARCH_ISA_A64)
template <>
constexpr KOKKOS_INLINE_FUNCTION bool
kk_is_a64fx_mem_space<Kokkos::HostSpace>() {
  return true;
}
#endif  // a64fx architectures

// Host function to determine free and total device memory.
// Will throw if execution space doesn't support this.
template <typename MemorySpace>
inline void kk_get_free_total_memory(size_t& /* free_mem */,
                                     size_t& /* total_mem */) {
  std::ostringstream oss;
  oss << "Error: memory space " << MemorySpace::name()
      << " does not support querying free/total memory.";
  throw std::runtime_error(oss.str());
}

#ifdef KOKKOS_ENABLE_CUDA
template <>
inline void kk_get_free_total_memory<Kokkos::CudaSpace>(size_t& free_mem,
                                                        size_t& total_mem) {
  cudaMemGetInfo(&free_mem, &total_mem);
}
template <>
inline void kk_get_free_total_memory<Kokkos::CudaUVMSpace>(size_t& free_mem,
                                                           size_t& total_mem) {
  cudaMemGetInfo(&free_mem, &total_mem);
}
template <>
inline void kk_get_free_total_memory<Kokkos::CudaHostPinnedSpace>(
    size_t& free_mem, size_t& total_mem) {
  cudaMemGetInfo(&free_mem, &total_mem);
}
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
inline void kk_get_free_total_memory<Kokkos::Experimental::HIPSpace>(
    size_t& free_mem, size_t& total_mem) {
  hipMemGetInfo(&free_mem, &total_mem);
}
#endif

// FIXME_SYCL Use compiler extension instead of low level interface when
// available. Also, we assume to query memory associated with the default queue.
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GPU)
template <>
inline void kk_get_free_total_memory<Kokkos::Experimental::SYCLDeviceUSMSpace>(
    size_t& free_mem, size_t& total_mem) {
  sycl::queue queue;
  sycl::device device = queue.get_device();
  auto level_zero_handle =
      sycl::get_native<sycl::backend::ext_oneapi_level_zero>(device);

  uint32_t n_memory_modules = 0;
  zesDeviceEnumMemoryModules(level_zero_handle, &n_memory_modules, nullptr);

  if (n_memory_modules != 1) {
    std::ostringstream oss;
    oss << "Error: number of memory modules for the SYCL backend: "
        << n_memory_modules
        << ". We only support querying free/total memory if exactly one memory "
           "module was found. Make sure that ZES_ENABLE_SYSMAN=1 is set at run "
           "time if no memeory modules were found!";
    throw std::runtime_error(oss.str());
  }

  zes_mem_handle_t memory_module_handle;
  zesDeviceEnumMemoryModules(level_zero_handle, &n_memory_modules,
                             &memory_module_handle);
  zes_mem_state_t memory_properties{
      ZES_STRUCTURE_TYPE_MEM_PROPERTIES,
  };
  zesMemoryGetState(memory_module_handle, &memory_properties);
  total_mem = memory_properties.size;
  free_mem  = memory_properties.free;
}

template <>
inline void kk_get_free_total_memory<Kokkos::Experimental::SYCLHostUSMSpace>(
    size_t& free_mem, size_t& total_mem) {
  kk_get_free_total_memory<Kokkos::Experimental::SYCLDeviceUSMSpace>(free_mem,
                                                                     total_mem);
}

template <>
inline void kk_get_free_total_memory<Kokkos::Experimental::SYCLSharedUSMSpace>(
    size_t& free_mem, size_t& total_mem) {
  kk_get_free_total_memory<Kokkos::Experimental::SYCLDeviceUSMSpace>(free_mem,
                                                                     total_mem);
}
#endif

template <typename ExecSpace>
inline int kk_get_max_vector_size() {
  return Kokkos::TeamPolicy<ExecSpace>::vector_length_max();
}

#ifdef KOKKOS_ENABLE_SYCL
template <>
inline int kk_get_max_vector_size<Kokkos::Experimental::SYCL>() {
  // FIXME SYCL: hardcoding to 8 is a workaround that seems to work for all
  // kernels. Wait for max subgroup size query to be fixed in SYCL and/or
  // Kokkos. Then TeamPolicy::vector_length_max() can be used for all
  // backends.
  return 8;
}
#endif

inline int kk_get_suggested_vector_size(const size_t nr, const size_t nnz,
                                        const ExecSpaceType exec_space) {
  int suggested_vector_size_ = 1;
  int max_vector_size        = 1;
  switch (exec_space) {
    case Exec_CUDA: max_vector_size = 32; break;
    case Exec_HIP: max_vector_size = 64; break;
    case Exec_SYCL:
      // FIXME SYCL: same as above - 8 is a workaround
      max_vector_size = 8;
      break;
    default:;
  }
  switch (exec_space) {
    default: break;
    case Exec_SERIAL:
    case Exec_OMP:
    case Exec_THREADS: break;
    case Exec_CUDA:
    case Exec_HIP:
    case Exec_SYCL:
      if (nr > 0) suggested_vector_size_ = nnz / double(nr) + 0.5;
      if (suggested_vector_size_ < 3) {
        suggested_vector_size_ = 2;
      } else if (suggested_vector_size_ <= 6) {
        suggested_vector_size_ = 4;
      } else if (suggested_vector_size_ <= 12) {
        suggested_vector_size_ = 8;
      } else if (suggested_vector_size_ <= 24) {
        suggested_vector_size_ = 16;
      } else if (suggested_vector_size_ <= 48) {
        suggested_vector_size_ = 32;
      } else {
        suggested_vector_size_ = 64;
      }
      if (suggested_vector_size_ > max_vector_size)
        suggested_vector_size_ = max_vector_size;
      break;
  }
  return suggested_vector_size_;
}

inline int kk_get_suggested_team_size(const int vector_size,
                                      const ExecSpaceType exec_space) {
  if (exec_space == Exec_CUDA || exec_space == Exec_HIP ||
      exec_space == Exec_SYCL) {
    // TODO: where this is used, tune the target value for
    // threads per block (but 256 is probably OK for CUDA and HIP)
    return 256 / vector_size;
  } else {
    return 1;
  }
}

// Taken from kokkos perf_test:
// https://github.com/kokkos/kokkos/blob/3e51447871eaec53ac4adc94d4d5376b7345b360/core/perf_test/PerfTest_ExecSpacePartitioning.cpp#L7-L36
namespace Experimental {

template <class ExecSpace>
struct SpaceInstance {
  static ExecSpace create() { return ExecSpace(); }
  static void destroy(ExecSpace&) {}
  static bool overlap() { return false; }
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct SpaceInstance<Kokkos::Cuda> {
  static Kokkos::Cuda create() {
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    return Kokkos::Cuda(stream);
  }
  static void destroy(Kokkos::Cuda& space) {
    cudaStream_t stream = space.cuda_stream();
    cudaStreamDestroy(stream);
  }
  static bool overlap() {
    bool value          = true;
    auto local_rank_str = std::getenv("CUDA_LAUNCH_BLOCKING");
    if (local_rank_str) {
      value = (std::atoi(local_rank_str) == 0);
    }
    return value;
  }
};
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
struct SpaceInstance<Kokkos::Experimental::HIP> {
  static Kokkos::Experimental::HIP create() {
    hipStream_t stream;
    hipStreamCreate(&stream);
    return Kokkos::Experimental::HIP(stream);
  }
  static void destroy(Kokkos::Experimental::HIP& space) {
    hipStream_t stream = space.hip_stream();
    hipStreamDestroy(stream);
  }
  static bool overlap() {
    // TODO: does HIP have an equivalent for CUDA_LAUNCH_BLOCKING?
    return true;
  }
};
#endif

}  // namespace Experimental

}  // namespace Impl
}  // namespace KokkosKernels

#endif
