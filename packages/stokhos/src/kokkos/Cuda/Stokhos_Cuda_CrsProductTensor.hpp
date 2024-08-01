// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_CUDA_CRS_PRODUCT_TENSOR_HPP

#include <iostream>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_CrsProductTensor.hpp"

#include "Stokhos_Cuda_DeviceProp.hpp"
#include "Stokhos_Cuda_WarpShuffle.hpp"

#include "Teuchos_TestForException.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

// Matrix-vector product specialization for CrsProductTensor layout
// To do:
//   * Incorporate texture/read-only cache
//   * Change tensor layout to allow coalesced reads with smaller
//     threads/row (will only probably help smaller problem sizes)
//   * Get average FEM entries/row from FEM graph (hardcoded to 27)
template< typename TensorScalar,
          typename MatrixScalar,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< TensorScalar, Kokkos::Cuda >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;

  typedef CrsProductTensor< TensorScalar, execution_space > tensor_type;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type;
  typedef Kokkos::View< VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda > vector_type;

#define USE_LDG 0

#if USE_LDG == 0

  // The multiply kernel
  class MultiplyKernel {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type BlockSize;

    MultiplyKernel( const matrix_type & A,
                    const vector_type & x,
                    const vector_type & y,
                    const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), BlockSize(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      // Get shared memory for loading x, A, and y
      volatile VectorScalar * const sh_x =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile MatrixScalar * const sh_A = sh_x + BlockSize*dim;
      volatile VectorScalar * const sh_y = sh_A + BlockSize*dim;
#if !HAVE_CUDA_SHUFFLE
      volatile VectorScalar * const sh_t = sh_y + dim;
#endif

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // Mask used for shuffle/warp-sync operations
      const int mask = blockDim.x == 32 ? 0xffffffff :
        ((1<<blockDim.x)-1)<<(threadIdx.y%(32/blockDim.x))*blockDim.x;

      // Zero y
      for ( size_type i = tid; i < dim; i += nid ) {
        sh_y[i] = 0.0;
      }

      // Loop over columns in the discrete (finite element) system.
      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      for (size_type iBlockEntry=iBlockEntryBeg; iBlockEntry<iBlockEntryEnd;
           iBlockEntry += BlockSize) {
        const size_type block_size =
          (iBlockEntryEnd-iBlockEntry < BlockSize) ?
            iBlockEntryEnd-iBlockEntry : BlockSize;

        // Wait for X and A to be used in the previous iteration
        // before reading new values.
        __syncthreads();

        // Coalesced read blocks of X and A into shared memory
        for ( size_type col = 0; col < block_size; ++col ) {

          const size_type iBlockColumn = m_A.graph.entries( iBlockEntry + col );
          const VectorScalar * const x = & m_x(        0, iBlockColumn );
          const MatrixScalar * const A = & m_A.values( 0, iBlockEntry + col );

          // Coalesced read by the whole block from global memory:
          for ( size_type i = tid; i < dim; i += nid ) {
            sh_x[col + i * BlockSize] = x[i]; // m_x(        i, iBlockColumn );
            sh_A[col + i * BlockSize] = A[i]; // m_A.values( i, iBlockEntry );
          }

        }

        __syncthreads(); // wait for X and A to be read before using them

        // This cuda block is responsible for computing all values of 'y'
        for ( size_type i = threadIdx.y; i < dim; i += blockDim.y ) {
          VectorScalar y = 0;

          // Product tensor entries which this warp will iterate:
          const size_type lBeg = m_A.block.entry_begin( i );
          const size_type lEnd = m_A.block.entry_end(   i );

          // Loop through sparse tensor contributions with coalesced reads.
          for ( size_type l = lBeg+threadIdx.x; l < lEnd; l += blockDim.x ) {

            // Read 'blockDim.x' entries from the tensor (coalesced)
            const size_type kj   = m_A.block.coord( l );
            const TensorScalar v = m_A.block.value( l );
            const size_type j    = ( kj & 0x0ffff ) * BlockSize ;
            const size_type k    = ( kj >> 16     ) * BlockSize ;

            for ( size_type col = 0; col < block_size; ++col ) {
              y += v * ( sh_A[col+j] * sh_x[col+k] +
                         sh_A[col+k] * sh_x[col+j] );
            }

          }

          // Reduction of 'y' within 'blockDim.x'
#if HAVE_CUDA_SHUFFLE
          if (blockDim.x >= 2) y += shfl_down(y, 1, blockDim.x, mask);
          if (blockDim.x >= 4) y += shfl_down(y, 2, blockDim.x, mask);
          if (blockDim.x >= 8) y += shfl_down(y, 4, blockDim.x, mask);
          if (blockDim.x >= 16) y += shfl_down(y, 8, blockDim.x, mask);
          if (blockDim.x >= 32) y += shfl_down(y, 16, blockDim.x, mask);
          if ( threadIdx.x == 0 ) sh_y[i] += y;
#else
          sh_t[ tid ] = y;
          if (threadIdx.x+16 < blockDim.x) sh_t[tid] += sh_t[tid+16];
          sync_warp(mask);
          if (threadIdx.x+ 8 < blockDim.x) sh_t[tid] += sh_t[tid+ 8];
          sync_warp(mask);
          if (threadIdx.x+ 4 < blockDim.x) sh_t[tid] += sh_t[tid+ 4];
          sync_warp(mask);
          if (threadIdx.x+ 2 < blockDim.x) sh_t[tid] += sh_t[tid+ 2];
          sync_warp(mask);
          if (threadIdx.x+ 1 < blockDim.x) sh_t[tid] += sh_t[tid+ 1];
          sync_warp(mask);
          if (threadIdx.x == 0) sh_y[i] += sh_t[tid];
#endif

        }

      }

      // Wait for all contributions of y to be completed
      __syncthreads();

      // Store result back in global memory
      for ( size_type i = tid; i < dim; i += nid ) {
        m_y( i, blockIdx.x ) = sh_y[ i ];
      }
    }
  };

#elif USE_LDG == 1

  // The multiply kernel -- using read-only data cache for A
  // Currently is slower than one above
  class MultiplyKernel {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type BlockSize;

    MultiplyKernel( const matrix_type & A,
                    const vector_type & x,
                    const vector_type & y,
                    const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), BlockSize(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      volatile VectorScalar * const sh_x =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile VectorScalar * const sh_y = sh_x + BlockSize*dim;
#if !HAVE_CUDA_SHUFFLE
      volatile VectorScalar * const sh_t = sh_y + dim;
#endif

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // Mask used for shuffle/warp-sync operations
      const int mask = blockDim.x == 32 ? 0xffffffff :
        ((1<<blockDim.x)-1)<<(threadIdx.y%(32/blockDim.x))*blockDim.x;

      // Zero y
      for ( size_type i = tid; i < dim; i += nid ) {
        sh_y[i] = 0.0;
      }

      // Loop over columns in the discrete (finite element) system.
      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      for (size_type iBlockEntry=iBlockEntryBeg; iBlockEntry<iBlockEntryEnd;
           iBlockEntry += BlockSize) {
        const size_type block_size =
          (iBlockEntryEnd-iBlockEntry < BlockSize) ?
            iBlockEntryEnd-iBlockEntry  : BlockSize;

        // Wait for X and A to be used in the previous iteration
        // before reading new values.
        __syncthreads();

        // Coalesced read blocks of X into shared memory
        for ( size_type col = 0; col < block_size; ++col ) {

          const size_type iBlockColumn = m_A.graph.entries( iBlockEntry + col );
          const VectorScalar * const x = & m_x( 0, iBlockColumn );

          // Coalesced read by the whole block from global memory:
          for ( size_type i = tid; i < dim; i += nid ) {
            sh_x[col + i * BlockSize] = x[i]; // m_x( i, iBlockColumn );
          }

        }

        __syncthreads(); // wait for X to be read before using them

        // This cuda block is responsible for computing all values of 'y'
        for ( size_type i = threadIdx.y; i < dim; i += blockDim.y ) {
          VectorScalar y = 0;

          // Product tensor entries which this warp will iterate:
          const size_type lBeg = m_A.block.entry_begin( i );
          const size_type lEnd = m_A.block.entry_end(   i );

          // Loop through sparse tensor contributions with coalesced reads.
          for ( size_type l = lBeg+threadIdx.x; l < lEnd; l += blockDim.x ) {

            // Read 'blockDim.x' entries from the tensor (coalesced)
            const size_type kj   = m_A.block.coord( l );
            const TensorScalar v = m_A.block.value( l );
            const size_type j    = ( kj & 0x0ffff ) ;
            const size_type k    = ( kj >> 16     ) ;

            for ( size_type col = 0; col < block_size; ++col ) {
              const size_type bCol = iBlockEntry + col;
#if (__CUDA_ARCH__ >= 350)
              y += v * ( __ldg(&m_A.values(j,bCol)) * sh_x[col+k*BlockSize] +
                         __ldg(&m_A.values(k,bCol)) * sh_x[col+j*BlockSize] );
#else
              y += v * ( m_A.values(j,bCol) * sh_x[col+k*BlockSize] +
                         m_A.values(k,bCol) * sh_x[col+j*BlockSize] );
#endif
            }

          }

          // Reduction of 'y' within 'blockDim.x'
#if HAVE_CUDA_SHUFFLE
          if (blockDim.x >= 2) y += shfl_down(y, 1, blockDim.x, mask);
          if (blockDim.x >= 4) y += shfl_down(y, 2, blockDim.x, mask);
          if (blockDim.x >= 8) y += shfl_down(y, 4, blockDim.x, mask);
          if (blockDim.x >= 16) y += shfl_down(y, 8, blockDim.x, mask);
          if (blockDim.x >= 32) y += shfl_down(y, 16, blockDim.x, mask);
          if ( threadIdx.x == 0 ) sh_y[i] += y;
#else
          sh_t[ tid ] = y;
          if (threadIdx.x+16 < blockDim.x) sh_t[tid] += sh_t[tid+16];
          sync_warp(mask);
          if (threadIdx.x+ 8 < blockDim.x) sh_t[tid] += sh_t[tid+ 8];
          sync_warp(mask);
          if (threadIdx.x+ 4 < blockDim.x) sh_t[tid] += sh_t[tid+ 4];
          sync_warp(mask);
          if (threadIdx.x+ 2 < blockDim.x) sh_t[tid] += sh_t[tid+ 2];
          sync_warp(mask);
          if (threadIdx.x+ 1 < blockDim.x) sh_t[tid] += sh_t[tid+ 1];
          sync_warp(mask);
          if (threadIdx.x == 0) sh_y[i] += sh_t[tid];
#endif

        }

      }

      // Wait for all contributions of y to be completed
      __syncthreads();

      // Store result back in global memory
      for ( size_type i = tid; i < dim; i += nid ) {
        m_y( i, blockIdx.x ) = sh_y[ i ];
      }
    }
  };

#endif

  //------------------------------------

  struct TensorReadEntry {
    size_type block_size, shmem, num_blocks, num_warp;
    double reads;
  };

  static void apply( const matrix_type & A,
                     const vector_type & x,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.extent(0) - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tensor_align = tensor_dimension;
    const size_type avg_tensor_entries_per_row = A.block.avg_entries_per_row();

    // Should compute this from FEM graph
    const size_type fem_nnz_per_row = 27;

    // Get device properties we need for whatever device is currently selected
    DeviceProp device_prop;
    const size_type shcap = device_prop.shared_memory_capacity;
    const size_type sh_granularity = device_prop.shared_memory_granularity;
    const size_type max_shmem_per_block = device_prop.max_shmem_per_block;
    const size_type max_blocks_per_sm = device_prop.max_blocks_per_sm;
    const size_type warp_size = device_prop.warp_size;
    const size_type warp_granularity = device_prop.warp_granularity;
    const size_type max_warps_per_block =
      std::min(device_prop.max_threads_per_block / warp_size,
               device_prop.max_warps_per_sm);
    const size_type min_warps_per_block = 1;
    const size_type max_regs_per_sm = device_prop.max_regs_per_sm;
    const size_type max_regs_per_block = device_prop.max_regs_per_block;
    const size_type reg_bank_size = device_prop.reg_bank_size;

    // Compute number of warps we can fit on each SM based on register limits
    // Use Cuda introspection to determine number of registers per thread
    //const size_type regs_per_thread = 46;
    const size_type regs_per_thread =
      device_prop.get_kernel_registers(
        Kokkos::Impl::cuda_parallel_launch_local_memory<MultiplyKernel>);
    const size_type regs_per_warp =
      (warp_size*regs_per_thread + reg_bank_size-1) & ~(reg_bank_size-1);
    const size_type warps_per_sm =
      (max_regs_per_sm/regs_per_warp) & ~(warp_granularity-1);
    const size_type warps_per_block =
      (max_regs_per_block/regs_per_warp) & ~(warp_granularity-1);

    // Compute number of threads per stochastic row based on number of
    // nonzero entries per row.
    // For double, 16 threads/row is still coalesced, but not for float.
    // We should reorder the tensor values for a given vector width to
    // maintain coalesced reads.  This would help smaller problems too by
    // allowing fewer threads/row.
    const size_type threads_per_row =
      avg_tensor_entries_per_row >= 88 ? 32 : 16;
    const size_type rows_per_warp = warp_size / threads_per_row;

    const size_type vec_scalar_size = sizeof(VectorScalar);
#if USE_LDG == 0
    const size_type mat_scalar_size = sizeof(MatrixScalar);
#endif

#define USE_FIXED_BLOCKSIZE 0

#if USE_FIXED_BLOCKSIZE

    const size_type num_blocks = 3;
    size_type nw = warps_per_sm / num_blocks;
    while (nw > 1 && num_blocks*nw % warp_granularity) --nw;
    const size_type num_warp = nw;
    const size_type sh_per_block = shcap / num_blocks;
    const size_type sr =
      device_prop.has_shuffle ? 0 : vec_scalar_size*warp_size*num_warp;
#if USE_LDG == 1
    size_type bs = ((sh_per_block - sr) / tensor_align - vec_scalar_size) /
      vec_scalar_size;
#else
    size_type bs = ((sh_per_block - sr) / tensor_align - vec_scalar_size) /
      (vec_scalar_size+mat_scalar_size);
#endif
    if (bs % 2 == 0) --bs;
    const size_type block_size_max = 31;
    const size_type block_size = std::min(bs, block_size_max);
    //const size_type block_size = 7;
#if USE_LDG == 1
    const size_type shmem =
      ( (vec_scalar_size*block_size+vec_scalar_size)*tensor_align + sr + sh_granularity-1 ) & ~(sh_granularity-1);
#else
    const size_type shmem =
      ( ((vec_scalar_size+mat_scalar_size)*block_size+vec_scalar_size)*tensor_align + sr + sh_granularity-1 ) & ~(sh_granularity-1);
#endif

#else

    // We want to maximize the number of blocks per SM (to maximize throughput)
    // as well as the block_size (to minimize tensor reads), subject to
    // shared memory and register constraints.  Here we try to do this computing
    // the number of tensor reads per block per thread for each choice of
    // block_size, and then choose the configuration with the smallest value.
    // This isn't perfect, but seems to generally work OK.  It could be
    // improved with a better model of:
    //   * Number of blocks versus warps per block (to minimize synchronization)
    //   * Thread efficiency for small numbers of rows per thread
    const size_type block_size_min = 3;
    const size_type half_nnz_per_row = fem_nnz_per_row / 2 + 1;
    const size_type block_size_max =
      half_nnz_per_row % 2 ? half_nnz_per_row + 1 : half_nnz_per_row;
    Teuchos::Array<TensorReadEntry> reads_per_thread;
    for (size_type bs = block_size_min; bs<=block_size_max; bs+=2) {
      // We don't know the number of warps yet, so we just have to bound
      // sr by the maximum number possible (which is all warps in 1 block)
      const size_type sr =
        device_prop.has_shuffle ? 0 : vec_scalar_size*warp_size*warps_per_block;
#if USE_LDG == 1
      size_type shmem =
        (vec_scalar_size*bs+vec_scalar_size)*tensor_align+sr;
#else
      size_type shmem =
        ((vec_scalar_size+mat_scalar_size)*bs+vec_scalar_size)*tensor_align+sr;
#endif
      shmem = (shmem + sh_granularity-1) & ~(sh_granularity-1);
      if (shmem <= max_shmem_per_block) {
        size_type num_blocks = std::min(shcap / shmem, max_blocks_per_sm);
        size_type tensor_reads = (fem_nnz_per_row+bs-1) / bs;
        size_type num_warp =
          std::min(std::max(std::min(warps_per_sm/num_blocks, warps_per_block),
                            min_warps_per_block),
                   max_warps_per_block);
        while (num_warp > 1 && num_blocks*num_warp % warp_granularity)
          --num_warp;
        TensorReadEntry entry;
        entry.block_size = bs;
        entry.shmem = shmem;
        entry.num_blocks = num_blocks;
        entry.num_warp = num_warp;

        // Prefer at least 3 blocks
        size_type factor = std::min(num_blocks,3u);
        entry.reads = (static_cast<double>(tensor_reads) /
                       static_cast<double>(factor*num_blocks*num_warp));
        reads_per_thread.push_back(entry);
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      reads_per_thread.size() == 0, std::logic_error,
      "Stochastic problem dimension is too large to fit in shared memory");
    size_type idx = 0;
    double min_reads = reads_per_thread[0].reads;
    for (int i=1; i<reads_per_thread.size(); ++i) {
      if (reads_per_thread[i].reads < min_reads) {
        idx = i;
        min_reads = reads_per_thread[i].reads;
      }
    }

    const size_type block_size = reads_per_thread[idx].block_size;
    const size_type shmem = reads_per_thread[idx].shmem;
    const size_type num_blocks = reads_per_thread[idx].num_blocks;
    const size_type num_warp = reads_per_thread[idx].num_warp;

#endif

    // Setup thread blocks and grid
    const dim3 dBlock( threads_per_row , rows_per_warp*num_warp , 1 );
    const dim3 dGrid( row_count, 1, 1 );

#if 0
    std::cout << "block_size = " << block_size
              << " tensor reads = " << (fem_nnz_per_row+block_size-1)/block_size
              << " regs_per_thread = " << regs_per_thread
              << " num blocks = " << num_blocks
              << " num warps = " << num_warp
              << " num rows = " << tensor_dimension
              << " rows/warp = " << tensor_dimension / (num_warp*rows_per_warp)
              << " avg entries/row = " <<  avg_tensor_entries_per_row
              << std::endl;
#endif

    // Finally launch our kernel
    //cudaProfilerStart();
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
      ( MultiplyKernel( A, x, y, block_size ) );
    //cudaProfilerStop();
  }

};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_CRS_PRODUCT_TENSOR_HPP */
