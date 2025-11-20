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

#if defined(STK_ENABLE_GPU)
#if !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
#define STK_ENABLE_GPU_BUT_NO_RDC
#endif
#endif

#define STK_INLINE_FUNCTION \
  STK_DEPRECATED_MSG("STK_INLINE_FUNCTION is deprecated. Use KOKKOS_INLINE_FUNCTION directly") KOKKOS_INLINE_FUNCTION
#define STK_FUNCTION STK_DEPRECATED_MSG("STK_FUNCTION is deprecated. Use KOKKOS_FUNCTION directly") KOKKOS_FUNCTION

// Until we can use the C++20 std::source_location capability, we must instead fall back
// on compiler extensions.  Detect what we have available.

#ifndef __has_builtin
  #define __has_builtin(x) 0  // Compatibility with non-LLVM-based compilers.
#endif

// GCC has the location built-ins
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
  #define HAS_GCC_LOCATION_BUILTINS
#endif

// Clang and the newer LLVM-based Intel compiler have the location build-ins
#if (defined(__clang__) || defined(__INTEL_LLVM_COMPILER)) && __has_builtin(__builtin_LINE)
  #define HAS_LLVM_LOCATION_BUILTINS
#endif

#if defined(HAS_GCC_LOCATION_BUILTINS) || defined(HAS_LLVM_LOCATION_BUILTINS)
  #define HAS_LOCATION_BUILTINS
#endif


#if defined(HAS_LOCATION_BUILTINS)
  #define STK_HOST_USE_LOCATION_BUILTINS
  #define STK_HOST_FILE __builtin_FILE()
  #define STK_HOST_LINE __builtin_LINE()
#else
  #define STK_HOST_FILE ""
  #define STK_HOST_LINE -1
#endif

// None of this is compatible with Clang or ROCm builds, so disable if we are doing a GPU-based device build
#if defined(HAS_LOCATION_BUILTINS) && !defined(STK_ENABLE_GPU)
  #define STK_DEVICE_USE_LOCATION_BUILTINS
  #define STK_DEVICE_FILE __builtin_FILE()
  #define STK_DEVICE_LINE __builtin_LINE()
#else
  #define STK_DEVICE_FILE ""
  #define STK_DEVICE_LINE -1
#endif

#endif
