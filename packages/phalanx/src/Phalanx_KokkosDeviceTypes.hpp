#ifndef PHALANX_KOKKOS_DEVICE_TYPES_HPP
#define PHALANX_KOKKOS_DEVICE_TYPES_HPP

//Kokkos includes
#include "Kokkos_Core.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Phalanx_config.hpp"

// ***************************************
// * DEVICE TYPE
// ***************************************

#if defined(PHX_KOKKOS_DEVICE_TYPE_CUDA)

//#include <Kokkos_Cuda.hpp>
namespace PHX { typedef Kokkos::Cuda Device; }
  
#elif defined(PHX_KOKKOS_DEVICE_TYPE_OPENMP)
  
//#include <Kokkos_hwloc.hpp>
//#include <Kokkos_OpenMP.hpp>
namespace PHX { typedef Kokkos::OpenMP Device; }

#elif defined(PHX_KOKKOS_DEVICE_TYPE_THREAD)
  
#include <Kokkos_hwloc.hpp>
//#include <Kokkos_Threads.hpp>
namespace PHX { typedef Kokkos::Threads Device; }

#elif defined(PHX_KOKKOS_DEVICE_TYPE_SERIAL)
  
//#include <Kokkos_Serial.hpp>
namespace PHX { typedef Kokkos::Serial Device; }

#endif

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

}

#endif
