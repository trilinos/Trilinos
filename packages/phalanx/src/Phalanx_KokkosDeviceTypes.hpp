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

  using exec_space = PHX::Device::execution_space;
  using mem_space  = PHX::Device::memory_space;

  using ExecSpace  = PHX::Device::execution_space;
  using MemSpace   = PHX::Device::memory_space;

}

// ***************************************
// * Kokkos View Properties
// ***************************************

namespace PHX {

  template <typename T> 
  struct remove_all_pointers{using type = T;};

  template <typename T> 
  struct remove_all_pointers<T*>{using type = typename PHX::remove_all_pointers<T>::type;};

  using DefaultDevLayout = PHX::exec_space::array_layout;

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
