#ifndef STK_KOKKOS_MACROS_H
#define STK_KOKKOS_MACROS_H

#include "Kokkos_Macros.hpp"

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#define STK_ENABLE_GPU
#endif

#if defined(STK_ENABLE_GPU)
  #ifndef STK_USE_DEVICE_MESH
    #define STK_USE_DEVICE_MESH
  #endif
#endif

#define STK_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#define STK_FUNCTION KOKKOS_FUNCTION

#endif

