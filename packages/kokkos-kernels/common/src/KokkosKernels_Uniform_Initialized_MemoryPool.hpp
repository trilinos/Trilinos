//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef _KOKKOSKERNELS_MEMPOOL_HPP
#define _KOKKOSKERNELS_MEMPOOL_HPP

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Utils.hpp"
#include <iostream>

namespace KokkosKernels {

namespace Impl {

enum PoolType { OneThread2OneChunk, ManyThread2OneChunk };

/*! \brief Class for simple memory pool allocations.
 *  At the constructor, we set the number of chunks and the size of a chunk,
 *  and afterwards each allocation returns a pointer in the size of 'chunk' that
 has
 *  been set in the constructor.
 *
 *  There are two modes of the memory pool:
 *  OneThread2OneChunk: This is the case where we have 1-to-1 chunks for the
 number of
 *    threads, that is #chunks = #threads. This is the usual case for CPUs,
 OpenMP and Threads.
 *    We have a dedicated memory for each thread.
 *
 *    NOTE: This mode can be used for GPUs as well, based on the maximum number
 of threads
 *    that can be run on GPUs. However, since we this mode require unique thread
 id's this
 *    is not yet can be retrieved on GPUs. Therefore, dont use this mode on GPUs
 yet.
 *
 *   ManyThread2OneChunk: This is the case where we have chunks <= # threads.
 Many thread
 *   will race for a memory allocation, some will get NULL pointers if there is
 not
 *   enough memory. This case still would work for #chunks = #threads, with an
 extra atomic
 *   operation. On GPUs, even when #chunks = Kokkos::Cuda().concurrency(), this
 option is safe
 *   to use.
 *
 *   This pool is limited in the sense that each allocation will be in the size
 of chunksize.
 *   The pool is specifically aimed for below cases which are difficult to write
 in kokkos:
 *
 *    #pragma omp parallel
 *    {
 *      ////////////initializations//////////
 *      for ( i: 1-> N) Accumulator[i] = 0;
 *      UsedAccumulatorIndices[N];
 *      UsedAccumulatorIndexCount = 0;
 *      #pragma omp for
 *      {
 *        //////work on accumulator//////
 *        ........
 *        ........
 *
 *        //////reset the accumulator//////
 *        for (i: 1->UsedAccumulatorIndexCount) Accumulator[
 UsedAccumulatorIndices[i] ] = 0;
 *      }
 *    }
 *
 *  In this example, each thread initializes its thread-private Accumulator
 (with size N) once.
 *  Then it clears its Accumulator at the end of each iteration using the sparse
 indices
 *  (which is usually much less than N).
 *  Doing that in Kokkos,
 *      --- either requires the initializations to go into loop body, resulting
 in N work in each loop.
 *      --- Or, we can use preinitialized 2d views where the first dimension is
 ExecutionSpace().concurrency()
 *          However, this case becomes a problem in CUDA, as concurrency is
 pretty high and we might not have
 *          enough memory for that.
 *
 *
 *
 *  Using this pool, on can simply use memory allocations without worrying about
 the memory constraints, with
 *  a sacrifize on the parallelization on cuda, when there is not enough memory.
 Also, for the algorithms that
 *  would use scratch space, and would need memory only if it runs out of
 scratch space, since this memory allocation
 *  will not be performed by all threads, there might not be sacrifize on the
 parallelism even when there is not
 *  enough memory for all threads.
 *
 *
 *  void operator()(const typename
 Kokkos::TeamPolicy<ExecutionSpace>::member_type & teamMember) const {
 *      volatile idx * myData = NULL;
 *      size_t tid = this->get_thread_id(teamMember);
 *
 *      while (myData == NULL){
 *        Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile idx *
 &memptr) {
 *          memptr = (volatile idx * )this->my_memory_pool.allocate_chunk(tid);
 *        }, myData);
 *      }
 *
 *
 *      /////..............work on memptr................../////
 *      /////..............reset as above on memptr so that thread leaves the
 memptr as it finds................../////
 *
 *      Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
          this->my_memory_pool.release_chunk((idx *) myData);
        });
 *  }
 *
 * It would only be set as below for different spaces.
 *  typedef typename KokkosKernels::Impl::UniformMemoryPool<ExecSpace,
 data_type> simple_pool;
 *  simple_pool mypool(number_of_chuns, chunksize, initialization_value,
 KokkosKernels::Impl::OneThread2OneChunk)
 *  or
 *  simple_pool mypool(number_of_chuns, chunksize, initialization_value,
 KokkosKernels::Impl::ManyThread2OneChunk)
 *
 *  Therefore, the option OneThread2OneChunk will guarantee that the same memory
 is used again and again by the same thread.
 *  This is not guaranteed by ManyThread2OneChunk, but as long as threads reset
 the memory they use, it will guarantee a memory
 *  that is has been initialized.
 */
template <typename MyExecSpace, typename data_type>
class UniformMemoryPool {
 private:
  typedef int lock_type;
  typedef typename Kokkos::View<lock_type *, MyExecSpace> lock_view_t;
  typedef typename Kokkos::View<data_type *, MyExecSpace> data_view_t;

  size_t num_chunks;
  size_t num_set_chunks;
  size_t modular_num_chunks;
  size_t chunk_size;
  size_t overall_size;
  // const size_t next_free_chunk;
  // const size_t last_free_chunk;

  // index_view_t free_chunks;
  lock_view_t chunk_locks;
  lock_type *pchunk_locks;
  data_view_t data_view;
  data_type *data;
  PoolType pool_type;

 public:
  using execution_space = typename MyExecSpace::execution_space;
  using memory_space    = typename MyExecSpace::memory_space;

  /**
   * \brief UniformMemoryPool constructor.
   * \param num_chunks_: number of chunks. This will be rounded to minimum pow2
   * number. \param chunk_size_: chunk size, the size of each allocation. \param
   * initialized_value: the value to initialize \param pool_type_: whether
   * ManyThread2OneChunk or OneThread2OneChunk
   */
  UniformMemoryPool(const size_t num_chunks_, const size_t set_chunk_size_, const data_type initialized_value = 0,
                    const PoolType pool_type_ = OneThread2OneChunk, bool initialize = true)
      : num_chunks(1),
        num_set_chunks(num_chunks_),
        modular_num_chunks(0),
        chunk_size(set_chunk_size_),
        overall_size(),
        // next_free_chunk(0),
        // last_free_chunk(chunk_size_),
        // free_chunks (),
        chunk_locks(),
        pchunk_locks(),
        data_view(),
        data(),
        pool_type(pool_type_) {
    num_chunks = 1;
    while (num_set_chunks > num_chunks) {
      num_chunks *= 2;
    }
    modular_num_chunks = num_chunks - 1;
    overall_size       = num_chunks * chunk_size;
    if (num_set_chunks > 0) {
      data_view = data_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "pool data"), overall_size);
    }
    data = (data_view.data());

    this->set_pool_type(pool_type_);

    if (initialize) {
      Kokkos::deep_copy(data_view, initialized_value);
    }
  }

  /**
   * \brief UniformMemoryPool constructor
   */
  UniformMemoryPool()
      : num_chunks(1),
        num_set_chunks(0),
        modular_num_chunks(0),
        chunk_size(0),
        overall_size(0),
        // next_free_chunk(0),
        // last_free_chunk(0),
        // free_chunks (),
        chunk_locks(),
        pchunk_locks(),
        data_view(),
        data(),
        pool_type() {}

  ~UniformMemoryPool() = default;

  UniformMemoryPool(UniformMemoryPool &&)                 = default;
  UniformMemoryPool(const UniformMemoryPool &)            = default;
  UniformMemoryPool &operator=(UniformMemoryPool &&)      = default;
  UniformMemoryPool &operator=(const UniformMemoryPool &) = default;

  /**
   * \brief To set the pool type
   * \param pool_type_: whether ManyThread2OneChunk or OneThread2OneChunk
   */
  void set_pool_type(PoolType pool_type_) {
    if (pool_type_ == ManyThread2OneChunk) {
      if (num_set_chunks) {
        chunk_locks = lock_view_t("locks", num_chunks);
      }
      pchunk_locks = chunk_locks.data();
    }
  }

  /**
   * \brief Print the content of memory pool
   */
  void print_memory_pool(bool print_all = false) const {
    std::cout << "num_chunks:" << num_chunks << std::endl;
    std::cout << "chunk_size:" << chunk_size << std::endl;
    std::cout << "overall_size:" << overall_size << std::endl;
    std::cout << "modular_num_chunks:" << modular_num_chunks << std::endl;

    // std::cout << "Printing free_chunks view" << std::endl;
    // print_1Dview(free_chunks, print_all);
    std::cout << "Printing chunk_locks view" << std::endl;
    print_1Dview(chunk_locks, print_all);
    std::cout << "Printing data view" << std::endl;
    print_1Dview(data_view, print_all);
  }

  /**
   * \brief Returns the unique memory location for thread.
   * This would uniquely return a memory location, for mode: OneThread2OneChunk
   * Use this function only if you dont wanna pay for the cost of switch in
   * allocate_chunk function, and you know that memory pool is always
   * OneThread2OneChunk. Otherwise you are likely to get a memory location that
   * is not private.
   */
  KOKKOS_INLINE_FUNCTION
  data_type *get_my_chunk(const size_t &thread_index) const {
    // return data + (thread_index % num_chunks) * chunk_size;
    return data + (thread_index & modular_num_chunks) * chunk_size;
  }

  /**
   * \brief Returns the unique memory location for thread.
   * This would uniquely return a memory location, for mode:
   * ManyThread2OneChunk. Use this function only if you dont wanna pay for the
   * cost of switch in allocate_chunk function, and you know that memory pool is
   * always ManyThread2OneChunk. Otherwise you will get a segfault.
   */
  KOKKOS_INLINE_FUNCTION
  data_type *get_arbitrary_free_chunk(const size_t &thread_index) const {
    return this->get_arbitrary_free_chunk(thread_index, num_chunks);
  }

  KOKKOS_INLINE_FUNCTION
  data_type *get_arbitrary_free_chunk(const size_t &thread_index, const size_t max_tries) const {
    size_t chunk_index = thread_index & modular_num_chunks;
    size_t num_try     = 0;
    while (!Kokkos::atomic_compare_exchange_strong(pchunk_locks + chunk_index, 0, 1)) {
      chunk_index = (chunk_index + 1) & modular_num_chunks;
      ++num_try;
      if (num_try > max_tries) {
        return NULL;
      }
    }
    return data + chunk_index * chunk_size;
  }

  /**
   * \brief Returns the unique memory location for thread.
   */
  KOKKOS_INLINE_FUNCTION
  data_type *allocate_chunk(const size_t &thread_index) const {
    switch (this->pool_type) {
      default:
      case OneThread2OneChunk:
        // printf("OneThread2OneChunk alloc for :%ld\n", thread_index);
        return this->get_my_chunk(thread_index);
      case ManyThread2OneChunk:
        // printf("ManyThread2OneChunk alloc for :%ld\n", thread_index);
        return this->get_arbitrary_free_chunk(thread_index, num_chunks);
    }
  }

  /**
   * \brief Releases the memory that has been allocated.
   * Use this function only if you dont wanna pay for the cost of switch in
   * release_chunk function, and you know that memory pool is always
   * ManyThread2OneChunk. Otherwise you will get a segfault.
   */
  KOKKOS_INLINE_FUNCTION
  void release_arbitrary_chunk(const data_type *chunk_ptr) const {
    size_t alloc_index = (chunk_ptr - data) / chunk_size;
    // printf("release:%ld #chunks:%ld\n", alloc_index, num_chunks);
    // chunk_locks(alloc_index) = false;
    chunk_locks(alloc_index) = 0;
  }

  /**
   * \brief Returns the chunk index of the pointer.
   */
  KOKKOS_INLINE_FUNCTION
  size_t get_chunk_index(const data_type *chunk_ptr) const { return (chunk_ptr - data) / chunk_size; }

  /**
   * \brief Releases the memory that has been allocated.
   */
  KOKKOS_INLINE_FUNCTION
  void release_chunk(const data_type *chunk_ptr) const {
    switch (this->pool_type) {
      default:
      case OneThread2OneChunk: break;
      case ManyThread2OneChunk: return this->release_arbitrary_chunk(chunk_ptr);
    }
  }
};

}  // namespace Impl
}  // namespace KokkosKernels

#endif
