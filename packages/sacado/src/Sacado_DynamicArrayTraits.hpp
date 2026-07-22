// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_DYNAMICARRAYTRAITS_HPP
#define SACADO_DYNAMICARRAYTRAITS_HPP

#include <new>
#include <cstring>
#include <stdint.h>

#include "Sacado_Traits.hpp"
#if defined(HAVE_SACADO_KOKKOS)
#include "Kokkos_Core.hpp"
#endif

namespace Sacado {

  template <typename ExecSpace>
  void createGlobalMemoryPool(const ExecSpace& space
            , const size_t min_total_alloc_size
            , const uint32_t min_block_alloc_size
            , const uint32_t max_block_alloc_size
            , const uint32_t min_superblock_size
            ) {}

  template <typename ExecSpace>
  void destroyGlobalMemoryPool(const ExecSpace& space) {}

#if 0 && defined(HAVE_SACADO_KOKKOS) && defined(KOKKOS_ENABLE_OPENMP)
  namespace Impl {
    extern const Kokkos::MemoryPool<Kokkos::OpenMP>* global_sacado_openmp_memory_pool;
  }

  inline void
  createGlobalMemoryPool(const ExecSpace& space
            , const size_t min_total_alloc_size
            , const uint32_t min_block_alloc_size
            , const uint32_t max_block_alloc_size
            , const uint32_t min_superblock_size
            )
  {
    typedef Kokkos::MemoryPool<Kokkos::OpenMP> pool_t;
    Impl::global_sacado_openmp_memory_pool =
      new pool_t(typename Kokkos::OpenMP::memory_space(),
          min_total_alloc_size,
          min_block_alloc_size,
          max_block_alloc_size,
          min_superblock_size);
  }

  inline void destroyGlobalMemoryPool(const Kokkos::OpenMP& space)
  {
    delete Impl::global_sacado_openmp_memory_pool;
  }
#endif

#if defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(KOKKOS_ENABLE_CUDA) && defined(__CUDACC__)

  namespace Impl {

    extern const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_host;
    extern const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_device;
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
    extern __device__ const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_on_device;
#else
    __device__ const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_on_device = 0;
#endif

    struct SetMemoryPoolPtr {
      Kokkos::MemoryPool<Kokkos::Cuda>* pool_device;
      __device__ inline void operator()(int) const {
        global_sacado_cuda_memory_pool_on_device = pool_device;
      };
    };

  }

  // For some reason we get memory errors if these functions are defined in
  // Sacado_DynamicArrayTraits.cpp
  inline void
  createGlobalMemoryPool(const Kokkos::Cuda& space
            , const size_t min_total_alloc_size
            , const uint32_t min_block_alloc_size
            , const uint32_t max_block_alloc_size
            , const uint32_t min_superblock_size
            )
  {
    typedef Kokkos::MemoryPool<Kokkos::Cuda> pool_t;
    pool_t* pool =
      new pool_t(typename Kokkos::Cuda::memory_space(),
          min_total_alloc_size,
          min_block_alloc_size,
          max_block_alloc_size,
          min_superblock_size);
    Impl::SetMemoryPoolPtr f;
    KOKKOS_IMPL_CUDA_SAFE_CALL( cudaMalloc( &f.pool_device, sizeof(pool_t) ) );
    KOKKOS_IMPL_CUDA_SAFE_CALL( cudaMemcpy( f.pool_device, pool,
                                sizeof(pool_t),
                                cudaMemcpyHostToDevice ) );
    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0,1),f);
    Impl::global_sacado_cuda_memory_pool_host = pool;
    Impl::global_sacado_cuda_memory_pool_device = f.pool_device;
  }

  inline void destroyGlobalMemoryPool(const Kokkos::Cuda& space)
  {
    KOKKOS_IMPL_CUDA_SAFE_CALL( cudaFree( (void*) Impl::global_sacado_cuda_memory_pool_device ) );
    delete Impl::global_sacado_cuda_memory_pool_host;
  }

#endif

#if !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(KOKKOS_ENABLE_CUDA) && defined(__CUDACC__)

  namespace Impl {

    // Compute warp lane/thread index
     __device__ inline int warpLane(const int warp_size = 32) {
      return ( threadIdx.x + (threadIdx.y + threadIdx.z*blockDim.y)*blockDim.x ) % warp_size;
    }

    // Reduce y across the warp and broadcast to all lanes
    template <typename T>
     __device__ inline T warpReduce(T y, const int warp_size = 32) {
      for (int i=1; i<warp_size; i*=2) {
        y += Kokkos::shfl_down(y, i, warp_size);
      }
      y = Kokkos::shfl(y, 0, warp_size);
      return y;
    }

    // Non-inclusive plus-scan up the warp, replacing the first entry with 0
    template <typename T>
    __device__ inline int warpScan(T y, const int warp_size = 32) {
      const int lane = warpLane();
      y = Kokkos::shfl_up(y, 1, warp_size);
      if (lane == 0)
        y = T(0);
      for (int i=1; i<warp_size; i*=2) {
        T t = Kokkos::shfl_up(y, i, warp_size);
        if (lane > i)
          y += t;
      }
      return y;
    }

    template <typename T>
    __device__ inline T warpBcast(T y, int id, const int warp_size = 32) {
      return Kokkos::shfl(y, id, warp_size);
    }

  }

#endif

  namespace Impl {

    template <typename T>
    SACADO_INLINE_FUNCTION
    static T* ds_alloc(const int sz) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
      T* m = 0;
      if (sz > 0)
        KOKKOS_IMPL_CUDA_SAFE_CALL( cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal ) );
#elif defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDA_ARCH__)
      // This code assumes all threads enter ds_alloc, even those with sz == 0
      T* m = 0;
      const int total_sz = warpReduce(sz);
      const int lane = warpLane();
      if (total_sz > 0 && lane == 0) {
        m = static_cast<T*>(global_sacado_cuda_memory_pool_on_device->allocate(total_sz*sizeof(T)));
        if (m == 0)
          Kokkos::abort("Allocation failed.  Kokkos memory pool is out of memory");
      }
      m = warpBcast(m,0);
      m += warpScan(sz);
#elif 0 && defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL) && defined(KOKKOS_ENABLE_OPENMP)
      T* m = 0;
      if (sz > 0) {
        if (global_sacado_openmp_memory_pool != 0) {
          m = static_cast<T*>(global_sacado_openmp_memory_pool->allocate(sz*sizeof(T)));
          if (m == 0)
            Kokkos::abort("Allocation failed.  Kokkos memory pool is out of memory");
        }
        else
          m = static_cast<T* >(operator new(sz*sizeof(T)));
      }
#else
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#if defined(HAVE_SACADO_KOKKOS)
        if (m == 0)
          Kokkos::abort("Allocation failed.");
#endif
      }
#endif
      return m;
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    static void ds_free(T* m, int sz) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
      if (sz > 0)
        KOKKOS_IMPL_CUDA_SAFE_CALL( cudaFree(m) );
#elif defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDA_ARCH__)
      const int total_sz = warpReduce(sz);
      const int lane = warpLane();
      if (total_sz > 0 && lane == 0) {
        global_sacado_cuda_memory_pool_on_device->deallocate((void*) m, total_sz*sizeof(T));
      }
#elif 0 && defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL) && defined(KOKKOS_ENABLE_OPENMP)
      if (sz > 0) {
        if (global_sacado_openmp_memory_pool != 0)
          global_sacado_openmp_memory_pool->deallocate((void*) m, sz*sizeof(T));
        else
          operator delete((void*) m);
      }
#else
      if (sz > 0)
        operator delete((void*) m);
#endif
    }

  }

  /*!
   * \brief Dynamic array allocation class that works for any type
   */
  template <typename T, bool isScalar = IsScalarType<T>::value>
  struct ds_array {

    //! Get memory for new array of length \c sz
    SACADO_INLINE_FUNCTION
    static T* get(int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      T* p = m;
      for (int i=0; i<sz; ++i)
        new (p++) T();
      return m;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      T* p = m;
      for (int i=0; i<sz; ++i)
        new (p++) T(0.0);
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(const T* src, int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      T* p = m;
      for (int i=0; i<sz; ++i)
        new (p++) T(*(src++));
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* strided_get_and_fill(const T* src, int stride, int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      T* p = m;
      for (int i=0; i<sz; ++i) {
        new (p++) T(*(src));
        src += stride;
      }
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void copy(const T* src, T*  dest, int sz) {
      for (int i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_copy(const T* src, int src_stride,
                                    T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = *(src);
        dest += dest_stride;
        src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      for (int i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = T(0.);
        dest += stride;
      }
    }

    //! Destroy array elements and release memory
    SACADO_INLINE_FUNCTION
    static void destroy_and_release(T* m, int sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
        b->~T();
      Impl::ds_free(m, sz);
    }
  };

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDA_ARCH__)

  namespace Impl {

    template <typename T>
    SACADO_INLINE_FUNCTION
    static T* ds_strided_alloc(const int sz) {
      T* m = 0;
      // Only do strided memory allocations when we are doing hierarchical
      // parallelism with a vector dimension of 32.  The limitation on the
      // memory pool allowing only a single thread in a warp to allocate
      // makes it too difficult to do otherwise.
      if (blockDim.x == 32) {
        //const int lane = warpLane();
        const int lane = threadIdx.x;
        if (sz > 0 && lane == 0) {
#if defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL)
          m = static_cast<T*>(global_sacado_cuda_memory_pool_on_device->allocate(sz*sizeof(T)));
          if (m == 0)
            Kokkos::abort("Allocation failed.  Kokkos memory pool is out of memory");
#else
          m = static_cast<T* >(operator new(sz*sizeof(T)));
#if defined(HAVE_SACADO_KOKKOS)
          if (m == 0)
            Kokkos::abort("Allocation failed.");
#endif
#endif
        }
        m = warpBcast(m,0,blockDim.x);
      }
      else {
        if (sz > 0) {
          m = static_cast<T* >(operator new(sz*sizeof(T)));
#if defined(HAVE_SACADO_KOKKOS)
          if (m == 0)
            Kokkos::abort("Allocation failed.");
#endif
        }
      }

      return m;
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    static void ds_strided_free(T* m, int sz) {
      if (blockDim.x == 32) {
        // const int lane = warpLane();
        const int lane = threadIdx.x;
        if (sz > 0 && lane == 0) {
#if defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL)
          global_sacado_cuda_memory_pool_on_device->deallocate((void*) m, sz*sizeof(T));
#else
          operator delete((void*) m);
#endif
        }
      }
      else {
        if (sz > 0)
          operator delete((void*) m);
      }

    }

  }

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T>
  struct ds_array<T,true> {

    //! Get memory for new array of length \c sz
    SACADO_INLINE_FUNCTION
    static T* get(int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      return m;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      for (int i=threadIdx.x; i<sz; i+=blockDim.x)
        m[i] = 0.0;
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(const T* src, int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      for (int i=threadIdx.x; i<sz; i+=blockDim.x)
        m[i] = src[i];
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* strided_get_and_fill(const T* src, int stride, int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      for (int i=threadIdx.x; i<sz; i+=blockDim.x)
        m[i] = src[i*stride];
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void copy(const T* src, T* dest, int sz) {
      for (int i=threadIdx.x; i<sz; i+=blockDim.x)
        dest[i] = src[i];
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_copy(const T* src, int src_stride,
                             T* dest, int dest_stride, int sz) {
      for (int i=threadIdx.x; i<sz; i+=blockDim.x) {
        dest[i*dest_stride] = src[i*src_stride];
      }
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      for (int i=threadIdx.x; i<sz; i+=blockDim.x)
        dest[i] = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_zero(T* dest, int stride, int sz) {
      for (int i=threadIdx.x; i<sz; i+=blockDim.x) {
        dest[i*stride] = T(0.);
      }
    }

    //! Destroy array elements and release memory
    SACADO_INLINE_FUNCTION
    static void destroy_and_release(T* m, int sz) {
      Impl::ds_strided_free(m, sz);
    }
  };

#elif defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDA_ARCH__)

  namespace Impl {

    template <typename T>
    SACADO_INLINE_FUNCTION
    static T* ds_strided_alloc(const int sz) {
      T* m = 0;
      // Only do strided memory allocations when we are doing hierarchical
      // parallelism with a vector dimension of 32.  The limitation on the
      // memory pool allowing only a single thread in a warp to allocate
      // makes it too difficult to do otherwise.
      if (blockDim.x == 32) {
        // const int total_sz = warpReduce(sz);
        // const int lane = warpLane();
        const int total_sz = warpReduce(sz, blockDim.x);
        const int lane = threadIdx.x;
        if (total_sz > 0 && lane == 0) {
#if defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL)
          m = static_cast<T*>(global_sacado_cuda_memory_pool_on_device->allocate(total_sz*sizeof(T)));
          if (m == 0)
            Kokkos::abort("Allocation failed.  Kokkos memory pool is out of memory");
#else
          m = static_cast<T* >(operator new(total_sz*sizeof(T)));
#if defined(HAVE_SACADO_KOKKOS)
          if (m == 0)
            Kokkos::abort("Allocation failed.");
#endif
#endif
        }
        m = warpBcast(m,0,blockDim.x);
        m += lane;
      }
      else {
        if (sz > 0) {
          m = static_cast<T* >(operator new(sz*sizeof(T)));
#if defined(HAVE_SACADO_KOKKOS)
          if (m == 0)
            Kokkos::abort("Allocation failed.");
#endif
        }
      }

      return m;
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    static void ds_strided_free(T* m, int sz) {
      if (blockDim.x == 32) {
        // const int total_sz = warpReduce(sz);
        // const int lane = warpLane();
        const int total_sz = warpReduce(sz, blockDim.x);
        const int lane = threadIdx.x;
        if (total_sz > 0 && lane == 0) {
#if defined(HAVE_SACADO_KOKKOS) && defined(SACADO_KOKKOS_USE_MEMORY_POOL)
          global_sacado_cuda_memory_pool_on_device->deallocate((void*) m, total_sz*sizeof(T));
#else
          operator delete((void*) m);
#endif
        }
      }
      else {
        if (sz > 0)
          operator delete((void*) m);
      }
    }
  }

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T>
  struct ds_array<T,true> {

    //! Get memory for new array of length \c sz
    SACADO_INLINE_FUNCTION
    static T* get(int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      return m;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      for (int i=0; i<sz; ++i)
        m[i*blockDim.x] = 0.0;
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(const T* src, int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      for (int i=0; i<sz; ++i)
        m[i*blockDim.x] = src[i*blockDim.x];
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* strided_get_and_fill(const T* src, int stride, int sz) {
      T* m = Impl::ds_strided_alloc<T>(sz);
      for (int i=0; i<sz; ++i)
        m[i*blockDim.x] = src[i*stride];
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void copy(const T* src, T* dest, int sz) {
      for (int i=0; i<sz; ++i)
        dest[i*blockDim.x] = src[i*blockDim.x];
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_copy(const T* src, int src_stride,
                             T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = *(src);
        dest += dest_stride;
        src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      for (int i=0; i<sz; ++i)
        dest[i*blockDim.x] = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = T(0.);
        dest += stride;
      }
    }

    //! Destroy array elements and release memory
    SACADO_INLINE_FUNCTION
    static void destroy_and_release(T* m, int sz) {
      Impl::ds_strided_free(m, sz);
    }
  };

#else

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T>
  struct ds_array<T,true> {

    //! Get memory for new array of length \c sz
    SACADO_INLINE_FUNCTION
    static T* get(int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      return m;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(int sz) {
      T* m = Impl::ds_alloc<T>(sz);
#if defined(__CUDACC__ ) || defined(__HIPCC__ )
      for (int i=0; i<sz; ++i)
        m[i] = 0.0;
#else
      if (sz > 0)
        std::memset(m,0,sz*sizeof(T));
#endif
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* get_and_fill(const T* src, int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      for (int i=0; i<sz; ++i)
        m[i] = src[i];
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    SACADO_INLINE_FUNCTION
    static T* strided_get_and_fill(const T* src, int stride, int sz) {
      T* m = Impl::ds_alloc<T>(sz);
      for (int i=0; i<sz; ++i)
        m[i] = src[i*stride];
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void copy(const T* src, T* dest, int sz) {
      if (sz > 0 && dest != NULL && src != NULL)
#if defined( __CUDACC__) || defined(__HIPCC__ )
        for (int i=0; i<sz; ++i)
          dest[i] = src[i];
#else
        std::memcpy(dest,src,sz*sizeof(T));
#endif
    }

    //! Copy array from \c src to \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_copy(const T* src, int src_stride,
                                    T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = *(src);
        dest += dest_stride;
        src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      if (sz > 0 && dest != NULL)
#if defined(__CUDACC__ ) || defined(__HIPCC__ )
        for (int i=0; i<sz; ++i)
          dest[i] = T(0.);
#else
        std::memset(dest,0,sz*sizeof(T));
#endif
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = T(0.);
        dest += stride;
      }
    }

    //! Destroy array elements and release memory
    SACADO_INLINE_FUNCTION
    static void destroy_and_release(T* m, int sz) {
      Impl::ds_free(m, sz);
    }
  };

#endif

} // namespace Sacado

#endif // SACADO_DYNAMICARRAY_HPP
