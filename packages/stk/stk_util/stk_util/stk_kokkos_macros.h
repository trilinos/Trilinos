#ifndef STK_KOKKOS_MACROS_H
#define STK_KOKKOS_MACROS_H

//When we are ready to depend on kokkos we will #include Kokkos_Macros.hpp
//and define our STK_FUNCTION macros in terms of the corresponding kokkos macros.
//Until then, do our own (duplicated) implementation.

#if defined(__CUDA_ARCH__)

#define STK_INLINE_FUNCTION __device__ __host__ inline
#define STK_FUNCTION __device__ __host__

#else

#define STK_INLINE_FUNCTION inline
#define STK_FUNCTION

#endif

#endif

