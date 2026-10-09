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

  static_assert(Kokkos::is_execution_space_v<PHX::DefaultExecutionSpace>,
                "Phalanx: the type configured through "
                "Phalanx_DEFAULT_EXECUTION_SPACE is not a Kokkos execution "
                "space.");
  static_assert(Kokkos::is_memory_space_v<PHX::DefaultMemorySpace>,
                "Phalanx: the memory space derived from "
                "Phalanx_DEFAULT_EXECUTION_SPACE is not a Kokkos memory "
                "space.");

  //! The execution space Phalanx runs in.
  using ExecutionSpace = PHX::Device::execution_space;

  //! The memory space Phalanx allocates in.
  using MemorySpace = PHX::Device::memory_space;

  // Shared space pairs a configured execution space with a memory space it
  // did not choose, so nothing but this stops a pairing that cannot work.
  // Kokkos would otherwise report it from somewhere deep inside the first
  // kernel launch.
  static_assert(Kokkos::SpaceAccessibility<PHX::ExecutionSpace,
                                           PHX::MemorySpace>::accessible,
                "Phalanx: the configured execution space cannot access the "
                "memory space it was paired with.  Check "
                "Phalanx_DEFAULT_EXECUTION_SPACE against "
                "Phalanx_ENABLE_SHARED_SPACE, which replaces the memory space "
                "with Kokkos::SharedSpace.");

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
  //
  // The stride follows PHX::ExecutionSpace, the space Phalanx actually runs
  // on, and NOT whichever backends happen to be compiled in.  A CUDA build
  // that sets Phalanx_DEFAULT_EXECUTION_SPACE to a host space has no vector
  // dimension and wants a stride of 1, but a chain of KOKKOS_ENABLE_ tests
  // would hand it 32 and silently disagree with the team policy.
  namespace Impl {
    template <typename T> inline constexpr bool dependent_false = false;

    template <typename ExecSpace>
    struct FadStride {
      static_assert(dependent_false<ExecSpace>,
                    "Phalanx: hierarchical parallelism is enabled but no FAD "
                    "stride is defined for this execution space.  The stride "
                    "must equal the vector_size passed to Kokkos::TeamPolicy, "
                    "so it cannot be guessed -- add a specialization above, or "
                    "build without Sacado_ENABLE_HIERARCHICAL / "
                    "Sacado_ENABLE_HIERARCHICAL_DFAD.");
    };

#if defined(KOKKOS_ENABLE_CUDA)
    template <> struct FadStride<Kokkos::Cuda> {static constexpr int value = 32;};
#endif
#if defined(KOKKOS_ENABLE_HIP)
    template <> struct FadStride<Kokkos::HIP> {static constexpr int value = 64;};
#endif
#if defined(KOKKOS_ENABLE_SYCL)
    template <> struct FadStride<Kokkos::SYCL> {static constexpr int value = 32;};
#endif
    // A host backend has no vector dimension to partition, so hierarchical is
    // a no-op there.  Carried anyway so the hierarchical code paths still
    // compile, including on a GPU build that runs Phalanx on the host.
#if defined(KOKKOS_ENABLE_SERIAL)
    template <> struct FadStride<Kokkos::Serial> {static constexpr int value = 1;};
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
    template <> struct FadStride<Kokkos::OpenMP> {static constexpr int value = 1;};
#endif
#if defined(KOKKOS_ENABLE_THREADS)
    template <> struct FadStride<Kokkos::Threads> {static constexpr int value = 1;};
#endif
  }

  using DefaultFadLayout =
    Sacado::LayoutContiguous<DefaultDevLayout,
                             PHX::Impl::FadStride<PHX::ExecutionSpace>::value>;

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
