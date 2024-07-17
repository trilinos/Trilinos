// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_Kokkos_Tools_CheckStreams.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include <limits>

// **********************************
// Ideally, we would like to also check allocations, deallocations,
// deep_copy, create_mirror, create_mirror_view, and
// create_mirror_view_and_copy as well. The kokkos tools will need to
// be modified to pass in the device id for these functions. For now
// we can only check parallel_* and fencing.
// **********************************

// Lambdas can only be converted to function pointers if they do not capture.
// Using a global non-static variable in an unnamed namespace to "capture" the
// device id.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
namespace {
  uint32_t phalanx_default_stream_device_id = std::numeric_limits<uint32_t>::max();

  void phalanx_kt_parallel_x_callback(char const *label, uint32_t device_id,
                                      uint64_t * /*kernel_id*/)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(device_id == phalanx_default_stream_device_id,
                               std::runtime_error,
                               "\"ERROR: the kernel \"" << label
                               << "\" with device id=" << device_id
                               << " is the same as the default stream id="
                             << phalanx_default_stream_device_id);
  }

  void phalanx_kt_fence_callback(char const *label, uint32_t device_id,
                                 uint64_t * /*fence_id*/)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(device_id == phalanx_default_stream_device_id,
                               std::runtime_error,
                               "\"ERROR: the fence \"" << label
                               << "\" with device id=" << device_id
                               << " is the same as the default stream id="
                               << phalanx_default_stream_device_id);
  }
}
#endif

void PHX::set_enforce_no_default_stream_use()
{
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  phalanx_default_stream_device_id = Kokkos::Tools::Experimental::device_id(PHX::Device());

  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(phalanx_kt_parallel_x_callback);
  Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(phalanx_kt_parallel_x_callback);
  Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(phalanx_kt_parallel_x_callback);
  Kokkos::Tools::Experimental::set_begin_fence_callback(phalanx_kt_fence_callback);
#endif
}

void PHX::unset_enforce_no_default_stream_use()
{
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(nullptr);
  Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(nullptr);
  Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(nullptr);
  Kokkos::Tools::Experimental::set_begin_fence_callback(nullptr);
#endif
}
