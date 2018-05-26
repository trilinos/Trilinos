#ifndef PHALANX_KOKKOS_DEVICE_TYPES_HPP
#define PHALANX_KOKKOS_DEVICE_TYPES_HPP

//Kokkos includes
#include "Kokkos_Core.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Phalanx_config.hpp"
#include "Sacado_Fad_ExpressionTraits.hpp"
#include <type_traits>

// ***************************************
// * DEVICE TYPE
// ***************************************

namespace PHX {

#if defined(PHX_KOKKOS_DEVICE_TYPE_CUDA)

  //#include <Kokkos_Cuda.hpp>
  typedef Kokkos::Cuda Device;

#elif defined(PHX_KOKKOS_DEVICE_TYPE_OPENMP)

  //#include <Kokkos_hwloc.hpp>
  //#include <Kokkos_OpenMP.hpp>
  typedef Kokkos::OpenMP Device;

#elif defined(PHX_KOKKOS_DEVICE_TYPE_THREAD)

#include <Kokkos_hwloc.hpp>
  //#include <Kokkos_Threads.hpp>
  typedef Kokkos::Threads Device;

#elif defined(PHX_KOKKOS_DEVICE_TYPE_SERIAL)

  //#include <Kokkos_Serial.hpp>
  typedef Kokkos::Serial Device;

#endif

  using exec_space = PHX::Device::execution_space;
  using mem_space = PHX::Device::memory_space;

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

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)

#if defined(KOKKOS_ENABLE_CUDA)
  // Contiguous layout with FAD stride of 32.  IMPORTANT: The
  // FadStride must be the same as the vector_size in the
  // Kokkos::TeamPolicy constructor. This value is only used for SFad
  // and SLFad, not for DFad.
  using DefaultFadLayout = Kokkos::LayoutContiguous<DefaultDevLayout,32>;
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
}

#endif
