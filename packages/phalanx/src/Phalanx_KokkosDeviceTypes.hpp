// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_KOKKOS_DEVICE_TYPES_HPP
#define PHALANX_KOKKOS_DEVICE_TYPES_HPP

//Kokkos includes
#include <Kokkos_Core_fwd.hpp>
#if !defined(KOKKOS_ENABLE_IMPL_VIEW_LEGACY)
#include "Sacado.hpp"
#else
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#endif
#include "Kokkos_Core.hpp"
#include "Phalanx_config.hpp"
#include "Sacado_Fad_ExpressionTraits.hpp"
#include <type_traits>

// ***************************************
// * DEVICE TYPE
// ***************************************

namespace PHX {

  // Phalanx_DEFAULT_EXECUTION_SPACE and Phalanx_DEFAULT_MEMORY_SPACE are
  // configured as free-form type names, so a plausible-looking mistake -- a
  // memory space in the execution space slot, say -- otherwise surfaces much
  // later as an unreadable template error.
  static_assert(Kokkos::is_execution_space_v<PHX::DefaultExecutionSpace>,
                "Phalanx: the type configured through "
                "Phalanx_DEFAULT_EXECUTION_SPACE is not a Kokkos execution "
                "space.");
  static_assert(Kokkos::is_memory_space_v<PHX::DefaultMemorySpace>,
                "Phalanx: the type configured through "
                "Phalanx_DEFAULT_MEMORY_SPACE is not a Kokkos memory space.");

  //! The execution space Phalanx runs in.
  using ExecutionSpace = PHX::Device::execution_space;

  //! The memory space Phalanx allocates in.
  using MemorySpace = PHX::Device::memory_space;

  // The two spaces are configured independently, so nothing but this stops a
  // pairing that cannot work -- a device memory space with a host execution
  // space, say.  Kokkos would otherwise report it from somewhere deep inside
  // the first kernel launch.
  static_assert(Kokkos::SpaceAccessibility<PHX::ExecutionSpace,
                                           PHX::MemorySpace>::accessible,
                "Phalanx: the configured execution space cannot access the "
                "configured memory space.  Check "
                "Phalanx_DEFAULT_EXECUTION_SPACE against "
                "Phalanx_DEFAULT_MEMORY_SPACE (or Phalanx_ENABLE_SHARED_SPACE, "
                "which sets the memory space to Kokkos::SharedSpace).");

  // Earlier names for the two types above.  Phalanx accumulated several
  // spellings of each; they all mean what ExecutionSpace and MemorySpace mean.
  // Deprecated -- see packages/phalanx/scripts/migrate_phx_device.py.
  // Configure with Phalanx_HIDE_DEPRECATED_CODE=ON to build as though they
  // were already gone, which is how to check that code no longer needs them.
#ifndef PHALANX_HIDE_DEPRECATED_CODE
  using exec_space PHALANX_DEPRECATED_MSG("Use PHX::ExecutionSpace") =
      PHX::ExecutionSpace;
  using mem_space PHALANX_DEPRECATED_MSG("Use PHX::MemorySpace") =
      PHX::MemorySpace;
  using ExecSpace PHALANX_DEPRECATED_MSG("Use PHX::ExecutionSpace") =
      PHX::ExecutionSpace;
  using MemSpace PHALANX_DEPRECATED_MSG("Use PHX::MemorySpace") =
      PHX::MemorySpace;
#endif

}

// ***************************************
// * Kokkos View Properties
// ***************************************

namespace PHX {

  template <typename T> 
  struct remove_all_pointers{using type = T;};

  template <typename T> 
  struct remove_all_pointers<T*>{using type = typename PHX::remove_all_pointers<T>::type;};

  using DefaultDevLayout = PHX::ExecutionSpace::array_layout;

#if defined(SACADO_GPU_HIERARCHICAL_DFAD) || defined(SACADO_GPU_HIERARCHICAL)

  // Contiguous layout whose FAD stride is the width of the vector
  // dimension: a warp on Cuda, a wavefront on HIP, a sub-group on SYCL.
  // IMPORTANT: The FadStride must be the same as the vector_size in the
  // Kokkos::TeamPolicy constructor. This value is only used for SFad and
  // SLFad, not for DFad.
#if defined(KOKKOS_ENABLE_CUDA)
  using DefaultFadLayout = Sacado::LayoutContiguous<DefaultDevLayout,32>;
#elif defined(KOKKOS_ENABLE_HIP)
  using DefaultFadLayout = Sacado::LayoutContiguous<DefaultDevLayout,64>;
#elif defined(KOKKOS_ENABLE_SYCL)
  using DefaultFadLayout = Sacado::LayoutContiguous<DefaultDevLayout,32>;
#elif defined(KOKKOS_ENABLE_SERIAL) || defined(KOKKOS_ENABLE_OPENMP) ||        \
      defined(KOKKOS_ENABLE_THREADS)
  // A host backend has no vector dimension to partition, so hierarchical is a
  // no-op here.  Carried anyway so the hierarchical code paths still compile
  // on a CPU-only build.
  using DefaultFadLayout = Sacado::LayoutContiguous<DefaultDevLayout,1>;
#else
#error "Phalanx: hierarchical parallelism is enabled but no FAD stride is defined for this backend.  The stride must equal the vector_size passed to Kokkos::TeamPolicy, so it cannot be guessed -- add a branch above for the new backend, or build without Sacado_ENABLE_HIERARCHICAL / Sacado_ENABLE_HIERARCHICAL_DFAD."
#endif

#else
  using DefaultFadLayout = DefaultDevLayout;
#endif

  template <typename DataType>
  struct DevLayout {
    using ScalarType = typename std::remove_const<typename PHX::remove_all_pointers<DataType>::type>::type;
    using type = typename std::conditional<Sacado::IsADType<ScalarType>::value,DefaultFadLayout,DefaultDevLayout>::type;
  };

  template<typename DataType>
  using View = Kokkos::View<DataType,typename PHX::DevLayout<DataType>::type,PHX::Device>;

  template<typename DataType>
  using AtomicView = Kokkos::View<DataType,typename PHX::DevLayout<DataType>::type,PHX::Device,Kokkos::MemoryTraits<Kokkos::Atomic>>;

  template<typename DataType>
  using UnmanagedView = Kokkos::View<DataType,typename PHX::DevLayout<DataType>::type,PHX::Device,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
}

// Hack for HIP compiler bug. Partial template specialization of class
// device functions incorrectly requires the __device__ flag.
#ifdef KOKKOS_ENABLE_HIP
#define PHALANX_HIP_HACK_KOKKOS_FUNCTION KOKKOS_FUNCTION
#else
#define PHALANX_HIP_HACK_KOKKOS_FUNCTION
#endif

#endif
