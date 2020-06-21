#ifndef STK_KOKKOS_MACROS_H
#define STK_KOKKOS_MACROS_H

#include <Kokkos_Macros.hpp>

// This should eventually need to be supplemented with checks for ROCM and
// other accelerator platforms
//
#ifdef KOKKOS_ENABLE_CUDA
  #ifndef STK_USE_DEVICE_MESH
    #define STK_USE_DEVICE_MESH
  #endif
#endif

#define STK_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#define STK_FUNCTION KOKKOS_FUNCTION

#endif

