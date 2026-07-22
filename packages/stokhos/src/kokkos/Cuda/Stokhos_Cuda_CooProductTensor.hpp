// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_COO_PRODUCT_TENSOR_HPP
#define STOKHOS_CUDA_COO_PRODUCT_TENSOR_HPP

#include <iostream>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_CooProductTensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

template< typename TensorScalar,
          typename MatrixScalar,
          typename VectorScalar,
          bool Pack>
class Multiply<
  BlockCrsMatrix< CooProductTensor< TensorScalar, Kokkos::Cuda, Pack >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;

  typedef CooProductTensor< TensorScalar, execution_space, Pack > tensor_type;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type;
  typedef Kokkos::View< VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda > vector_type;

  typedef int rows_type;
  static const rows_type invalid_row = -1;

  class CooKernel {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type m_block_size;
    const size_type m_entries_per_thread;

    CooKernel( const matrix_type & A,
               const vector_type & x,
               const vector_type & y,
               const size_type entries_per_thread,
               const size_type block_size )
      : m_A(A),
        m_x(x),
        m_y(y),
        m_entries_per_thread(entries_per_thread),
        m_block_size(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type num_entries = m_A.block.entry_count();

      // Thread dimensions and index
      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // Shared memory
      VectorScalar * const sh_x =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_A = sh_x + m_block_size*dim;
      VectorScalar * const sh_y = sh_A + m_block_size*dim;
      volatile VectorScalar * const vals = sh_y + dim;
      volatile rows_type * const rows =
        reinterpret_cast<volatile rows_type*>(vals + nid);

      // Product tensor entries which this warp will iterate:
      const size_type entries_per_warp = blockDim.x * m_entries_per_thread;
      const size_type lBeg = threadIdx.y * entries_per_warp + threadIdx.x;
      const size_type lEnd = ( lBeg + entries_per_warp < num_entries ?
                               lBeg + entries_per_warp : num_entries );

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type femBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type femEnd = m_A.graph.row_map[ blockIdx.x + 1 ];

      // Zero y
      for (size_type l = tid; l < dim; l += nid) {
        sh_y[l] = 0.0;
      }

      // Initialize rows & vals arrays
      rows[tid] = invalid_row;
      vals[tid] = 0.0;

      // Loop over columns in the discrete (finite element) system.
      for (size_type femEntry=femBeg; femEntry<femEnd; femEntry+=m_block_size) {
        const size_type block_size =
          femEntry + m_block_size < femEnd ? m_block_size : femEnd - femEntry;

        // Wait for X and A to be used in the previous iteration
        // before reading new values.
        __syncthreads();

        // Coalesced read blocks of X and A into shared memory
        for (size_type col = 0; col < block_size; ++col) {

          const size_type femColumn = m_A.graph.entries( femEntry + col );
          const VectorScalar * const x = & m_x( 0, femColumn );
          const MatrixScalar * const A = & m_A.values( 0, femEntry + col );

          // Coalesced read by the whole block from global memory:
          for (size_type l = tid; l < dim; l += nid) {
            sh_x[col + l * m_block_size] = x[l];
            sh_A[col + l * m_block_size] = A[l];
          }

        }

        __syncthreads(); // wait for X and A to be read before using them

        // This cuda block is responsible for computing all values of 'y'
        for (size_type l = lBeg; l < lEnd; l += blockDim.x) {

          // Read 'blockDim.x' entries from the tensor (coalesced)
          size_type i, j, k;
          m_A.block.coord(l,i,j,k);
          const TensorScalar v = m_A.block.value( l );
          j *= m_block_size;
          k *= m_block_size;

          // Register storing local accumulation for row i
          VectorScalar y = 0;

          // Check for carry from previous row
          if (threadIdx.x == 0) {
            if  (i == rows[tid+31])
              y += vals[tid+31];
            else
              sh_y[rows[tid+31]] += vals[tid+31];
          }

          // Accumulate local row for the set of FEM columns
          for (size_type col = 0; col < block_size; ++col) {
            y += v * ( sh_A[col+j] * sh_x[col+k] + sh_A[col+k] * sh_x[col+j] );
          }

          // Store row and value into shared arrays
          rows[tid] = i;
          vals[tid] = y;

          // Reduce 'y' within 'blockDim.x' to the right for threads
          // on the same row
          if (threadIdx.x >=  1 && i == rows[tid- 1]) vals[tid] += vals[tid- 1];
          if (threadIdx.x >=  2 && i == rows[tid- 2]) vals[tid] += vals[tid- 2];
          if (threadIdx.x >=  4 && i == rows[tid- 4]) vals[tid] += vals[tid- 4];
          if (threadIdx.x >=  8 && i == rows[tid- 8]) vals[tid] += vals[tid- 8];
          if (threadIdx.x >= 16 && i == rows[tid-16]) vals[tid] += vals[tid-16];

          // Add local accumulation of y into sh_y for threads on
          // distinct rows.  We don't store thread 31 and instead carry it
          // to the next iteration to eliminate race conditions between warps
          if (threadIdx.x < 31 && i != rows[tid + 1])
            sh_y[i] += vals[tid];

        }

        // At this point we have blockDim.y values that need to be
        // reduced and stored.  Move these row/value pairs to the beginning
        __syncthreads();
        if (threadIdx.x == 31) {
          rows[threadIdx.y] = rows[tid];
          vals[threadIdx.y] = vals[tid];
        }
        __syncthreads();

        // Reduce these values to the right using the first warp
        // This assumes blockDim.x >= blockDim.y
        if (threadIdx.y == 0 && threadIdx.x < blockDim.y) {
          const size_type i = rows[tid];
          if (threadIdx.x >=  1 && i == rows[tid- 1]) vals[tid] += vals[tid- 1];
          if (threadIdx.x >=  2 && i == rows[tid- 2]) vals[tid] += vals[tid- 2];
          if (threadIdx.x >=  4 && i == rows[tid- 4]) vals[tid] += vals[tid- 4];
          if (threadIdx.x >=  8 && i == rows[tid- 8]) vals[tid] += vals[tid- 8];
          if (threadIdx.x >= 16 && i == rows[tid-16]) vals[tid] += vals[tid-16];

          if ((threadIdx.x == blockDim.y-1) ||
              (threadIdx.x  < blockDim.y-1 && i != rows[tid+1]))
            sh_y[i] += vals[tid];
        }

        // Reset rows and vals to prohibit carry across FEM columns
        rows[tid] = invalid_row;
        vals[tid] = 0.0;

      }

      // Wait for all contributions of y to be completed
      __syncthreads();

      // Store result back in global memory
      for (size_type l = tid; l < dim; l += nid) {
        m_y( l, blockIdx.x ) = sh_y[ l ];
      }
    }
  };

  //------------------------------------

  static void apply( const matrix_type & A,
                     const vector_type & x,
                     const vector_type & y )
  {
    const size_type fem_rows = A.graph.row_map.extent(0) - 1;
    const size_type stoch_rows = A.block.dimension();
    const size_type stoch_entries = A.block.entry_count();
    const size_type warp_size = Kokkos::Impl::CudaTraits::WarpSize;

#ifdef STOKHOS_DEBUG
    const size_type num_warp_max = 16; // Use fewer warps in debug mode to prevent
                                       // launch failures
#else
    const size_type num_warp_max = 20;
#endif
    const size_type num_warp =
      std::min( num_warp_max, stoch_entries / warp_size );
    const dim3 dBlock( warp_size , num_warp, 1 );
    const dim3 dGrid( fem_rows, 1, 1 );

    const size_type num_thread = dBlock.x*dBlock.y;
    const size_type entries_per_thread =
      (stoch_entries + num_thread-1) / num_thread;

    // Use at most half of shared memory to get 2 blocks per SMP
    const size_type size_rows = sizeof(rows_type) * num_thread;
    const size_type size_vals = sizeof(VectorScalar) * num_thread;
    const size_type shcap =
      Kokkos::Cuda().impl_internal_space_instance()->m_deviceProp.sharedMemPerBlock / 2;
    size_type bs =
      ((shcap-size_rows-size_vals) / (sizeof(VectorScalar)*stoch_rows) - 1) / 2;
    if (bs % 2 == 0) --bs; // Make block-size odd to reduce bank conflicts
    const size_type block_size_max = 31;
    const size_type block_size = std::min(bs, block_size_max);
    const size_type shmem =
      sizeof(VectorScalar) * ((2*block_size+1) * stoch_rows) + // A, x, y
      size_vals + size_rows;

#if 0
    //std::cout << std::endl << A.block << std::endl;
    const size_type fem_nnz = A.values.extent(1);
    std::cout << "Multiply< BlockCrsMatrix< CooProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  fem_rows(" << fem_rows << ")" << std::endl
              << "  fem_nnz(" << fem_nnz << ")" << std::endl
              << "  stoch_rows(" << stoch_rows << ")" << std::endl
              << "  stoch_nnz(" << stoch_entries << ")" << std::endl
              << "  threads_per_block(" << num_thread << ")" << std::endl
              << "  entries_per_thread(" << entries_per_thread << ")" << std::endl
             ;
#endif

    //cudaProfilerStart();
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
      ( CooKernel( A, x, y, entries_per_thread, block_size ) );
    //cudaProfilerStop();
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_COO_PRODUCT_TENSOR_HPP */
