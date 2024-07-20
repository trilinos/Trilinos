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
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Kokkos_Core.hpp"
#include "Phalanx_config.hpp"
#include "Sacado_Fad_ExpressionTraits.hpp"
#include <type_traits>

// ***************************************
// * DEVICE TYPE
// ***************************************

namespace PHX {

#if defined(PHX_KOKKOS_DEVICE_TYPE_CUDA)
  using Device = Kokkos::Cuda;
#elif defined(PHX_KOKKOS_DEVICE_TYPE_HIP)
  using Device = Kokkos::HIP;
#elif defined(PHX_KOKKOS_DEVICE_TYPE_SYCL)
  using Device = Kokkos::Experimental::SYCL;
#elif defined(PHX_KOKKOS_DEVICE_TYPE_OPENMP)
  using Device = Kokkos::OpenMP;
#elif defined(PHX_KOKKOS_DEVICE_TYPE_THREAD)
  #include <Kokkos_hwloc.hpp>
  using Device = Kokkos::Threads;
#elif defined(PHX_KOKKOS_DEVICE_TYPE_SERIAL)
  using Device = Kokkos::Serial;
#endif

  using exec_space = PHX::Device::execution_space;
  using mem_space  = PHX::Device::memory_space;

  using ExecSpace  = PHX::Device::execution_space;
  using MemSpace   = PHX::Device::memory_space;

}

// ***************************************
// * INDEX SIZE TYPE
// ***************************************

namespace PHX {

#if defined(PHX_INDEX_SIZE_TYPE_KOKKOS)
  typedef PHX::Device::size_type index_size_type;
#elif defined(PHX_INDEX_SIZE_TYPE_INT)
  typedef int index_size_type;
#elif defined(PHX_INDEX_SIZE_TYPE_UINT)
  typedef unsigned int index_size_type;
#elif defined(PHX_INDEX_SIZE_TYPE_LONGINT)
  typedef long int index_size_type;
#elif defined(PHX_INDEX_SIZE_TYPE_ULONGINT)
  typedef unsigned long int index_size_type;
#endif

  using index_t = index_size_type;

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

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD) || defined(SACADO_VIEW_CUDA_HIERARCHICAL)

  // Contiguous layout with FAD stride of 32 for cuda warp of 64 for
  // HIP warp.  IMPORTANT: The FadStride must be the same as the
  // vector_size in the Kokkos::TeamPolicy constructor. This value is
  // only used for SFad and SLFad, not for DFad.
#if defined(KOKKOS_ENABLE_CUDA)
  using DefaultFadLayout = Kokkos::LayoutContiguous<DefaultDevLayout,32>;
#elif defined(KOKKOS_ENABLE_HIP)
  using DefaultFadLayout = Kokkos::LayoutContiguous<DefaultDevLayout,64>;
#else
  using DefaultFadLayout = Kokkos::LayoutContiguous<DefaultDevLayout,1>;
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
