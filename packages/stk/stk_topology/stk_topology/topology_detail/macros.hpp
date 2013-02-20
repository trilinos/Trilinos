#ifndef STKTOPOLOGY_DETAIL_MACROS_HPP
#define STKTOPOLOGY_DETAIL_MACROS_HPP

#if defined( __CUDACC__ ) && defined( __CUDA_ARCH__ )

#define STKTOPOLOGY_INLINE_FUNCTION  __device__  __host__  inline

#else

#define STKTOPOLOGY_INLINE_FUNCTION inline

#endif

#endif //STKTOPOLOGY_DETAIL_MACROS_HPP
