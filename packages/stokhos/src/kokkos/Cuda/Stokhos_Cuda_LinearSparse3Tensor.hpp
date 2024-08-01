// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_LINEAR_SPARSE_3_TENSOR_HPP
#define STOKHOS_CUDA_LINEAR_SPARSE_3_TENSOR_HPP

#include <iostream>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_LinearSparse3Tensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

/*
 * This version uses 4 warps per block, one warp per spatial row, and loops
 * through stochastic rows 32 at a time.
 */
template< typename TensorScalar,
          typename MatrixScalar,
          typename VectorScalar,
          int BlockSize >
class Multiply<
  BlockCrsMatrix< LinearSparse3Tensor< TensorScalar, Kokkos::Cuda, BlockSize >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;

  typedef LinearSparse3Tensor< TensorScalar, execution_space, BlockSize > tensor_type;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type;
  typedef Kokkos::View< VectorScalar**,
                             Kokkos::LayoutLeft,
                             Kokkos::Cuda > vector_type;

  template <int MAX_COL>
  class ApplyKernelSymmetric {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;

    ApplyKernelSymmetric( const matrix_type & A,
                          const vector_type & x,
                          const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      volatile VectorScalar * const sh =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile VectorScalar * const sh_y0 =
        sh + blockDim.x*threadIdx.y;
      volatile VectorScalar * const sh_a0 =
        sh + blockDim.x*blockDim.y + MAX_COL*threadIdx.y;
      volatile VectorScalar * const sh_x0 =
        sh + blockDim.x*blockDim.y + MAX_COL*blockDim.y + MAX_COL*threadIdx.y;
      volatile size_type * const sh_col =
        (size_type*)(sh + blockDim.x*blockDim.y + 2*MAX_COL*blockDim.y) + MAX_COL*threadIdx.y;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type row = blockIdx.x*blockDim.y + threadIdx.y;
      if (row < m_A.graph.row_map.extent(0)-1) {
        const size_type colBeg = m_A.graph.row_map[ row ];
        const size_type colEnd = m_A.graph.row_map[ row + 1 ];

        // Read tensor entries
        const TensorScalar c0 = m_A.block.value(0);
        const TensorScalar c1 = m_A.block.value(1);

        // Zero y
        VectorScalar y0 = 0.0;

        // Read column offsets for this row
        for ( size_type lcol = threadIdx.x; lcol < colEnd-colBeg;
              lcol += blockDim.x )
          sh_col[lcol] = m_A.graph.entries( lcol+colBeg );

        // Loop over stochastic rows
        for ( size_type stoch_row = threadIdx.x; stoch_row < dim;
              stoch_row += blockDim.x) {

          VectorScalar yi = 0.0;

          // Loop over columns in the discrete (finite element) system.
          // We should probably only loop over a portion at a time to
          // make better use of L2 cache.
          //
          // We could also apply some warps to columns in parallel.
          // This requires inter-warp reduction of y0 and yi.
          for ( size_type col_offset = colBeg; col_offset < colEnd;
                col_offset += 1 ) {

            // Get column index
            const size_type lcol = col_offset-colBeg;
            const size_type col = sh_col[lcol];

            // Read a[i] and x[i] (coalesced)
            const MatrixScalar ai = m_A.values( stoch_row, col_offset );
            const VectorScalar xi = m_x( stoch_row, col );

            // Store a[0] and x[0] to shared memory for this column
            if (stoch_row == 0) {
              sh_a0[lcol] = ai;
              sh_x0[lcol] = xi;
            }

            // Retrieve a[0] and x[0] from shared memory for this column
            const MatrixScalar a0 = sh_a0[lcol];
            const VectorScalar x0 = sh_x0[lcol];

            // Subtract off contribution of first iteration of loop
            if (stoch_row == 0) y0 += (c0-3.0*c1)*a0*x0;

            yi += c1*(a0*xi + ai*x0);
            y0 += c1*ai*xi;
          }

          // Write y back to global memory
          m_y( stoch_row, row ) = yi;

          // Sync warps so rows don't get too out-of-sync to make use of
          // locality in spatial matrix and L2 cache
          __syncthreads();
        }

        // Reduction of 'y0' within warp
        sh_y0[ threadIdx.x ] = y0;
        if ( threadIdx.x + 16 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+16];
        if ( threadIdx.x +  8 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 8];
        if ( threadIdx.x +  4 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 4];
        if ( threadIdx.x +  2 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 2];
        if ( threadIdx.x +  1 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 1];

        // Store y0 back in global memory
        if ( threadIdx.x == 0 ) m_y( 0, row ) += sh_y0[0];

      }
    }
  };

  template <int MAX_COL>
  class ApplyKernelAsymmetric {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;

    ApplyKernelAsymmetric( const matrix_type & A,
                           const vector_type & x,
                           const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      volatile VectorScalar * const sh =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile VectorScalar * const sh_y0 =
        sh + blockDim.x*threadIdx.y;
      volatile VectorScalar * const sh_a0 =
        sh + blockDim.x*blockDim.y + MAX_COL*threadIdx.y;
      volatile VectorScalar * const sh_x0 =
        sh + blockDim.x*blockDim.y + MAX_COL*blockDim.y + MAX_COL*threadIdx.y;
      volatile size_type * const sh_col =
        (size_type*)(sh + blockDim.x*blockDim.y + 2*MAX_COL*blockDim.y) + MAX_COL*threadIdx.y;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type row = blockIdx.x*blockDim.y + threadIdx.y;
      if (row < m_A.graph.row_map.extent(0)-1) {
        const size_type colBeg = m_A.graph.row_map[ row ];
        const size_type colEnd = m_A.graph.row_map[ row + 1 ];

        // Read tensor entries
        const TensorScalar c0 = m_A.block.value(0);
        const TensorScalar c1 = m_A.block.value(1);
        const TensorScalar c2 = m_A.block.value(2);

        // Zero y
        VectorScalar y0 = 0.0;

        // Read column offsets for this row
        for ( size_type lcol = threadIdx.x; lcol < colEnd-colBeg;
              lcol += blockDim.x )
          sh_col[lcol] = m_A.graph.entries( lcol+colBeg );

        // Loop over stochastic rows
        for ( size_type stoch_row = threadIdx.x; stoch_row < dim;
              stoch_row += blockDim.x) {

          VectorScalar yi = 0.0;

          // Loop over columns in the discrete (finite element) system.
          // We should probably only loop over a portion at a time to
          // make better use of L2 cache.
          //
          // We could also apply some warps to columns in parallel.
          // This requires inter-warp reduction of y0 and yi.
          for ( size_type col_offset = colBeg; col_offset < colEnd;
                col_offset += 1 ) {

            // Get column index
            const size_type lcol = col_offset-colBeg;
            const size_type col = sh_col[lcol];

            // Read a[i] and x[i] (coalesced)
            const MatrixScalar ai = m_A.values( stoch_row, col_offset );
            const VectorScalar xi = m_x( stoch_row, col );

            // Store a[0] and x[0] to shared memory for this column
            if (stoch_row == 0) {
              sh_a0[lcol] = ai;
              sh_x0[lcol] = xi;
            }

            // Retrieve a[0] and x[0] from shared memory for this column
            const MatrixScalar a0 = sh_a0[lcol];
            const VectorScalar x0 = sh_x0[lcol];

            // Subtract off contribution of first iteration of loop
            if (stoch_row == 0) y0 += (c0-3.0*c1-c2)*a0*x0;

            yi += c1*(a0*xi + ai*x0) + c2*ai*xi;
            y0 += c1*ai*xi;
          }

          // Write y back to global memory
          m_y( stoch_row, row ) = yi;

          // Sync warps so rows don't get too out-of-sync to make use of
          // locality in spatial matrix and L2 cache
          __syncthreads();
        }

        // Reduction of 'y0' within warp
        sh_y0[ threadIdx.x ] = y0;
        if ( threadIdx.x + 16 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+16];
        if ( threadIdx.x +  8 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 8];
        if ( threadIdx.x +  4 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 4];
        if ( threadIdx.x +  2 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 2];
        if ( threadIdx.x +  1 < blockDim.x )
          sh_y0[threadIdx.x] += sh_y0[threadIdx.x+ 1];

        // Store y0 back in global memory
        if ( threadIdx.x == 0 ) m_y( 0, row ) += sh_y0[0];

      }
    }
  };


  //------------------------------------

  static void apply( const matrix_type & A,
                     const vector_type & x,
                     const vector_type & y )
  {
    const size_type num_spatial_rows = A.graph.row_map.extent(0) - 1;
    const size_type num_stoch_rows = A.block.dimension();

    const size_type warp_size = Kokkos::Impl::CudaTraits::WarpSize;
    const size_type num_rows_per_block = 4;
    size_type num_blocks = num_spatial_rows / num_rows_per_block;
    if (num_blocks * num_rows_per_block < num_spatial_rows)
      ++num_blocks;
    const dim3 dBlock( warp_size, num_rows_per_block );
    const dim3 dGrid( num_blocks , 1, 1 );

    const int MAX_COL = 32; // Maximum number of nonzero columns per row
    const size_type shmem =
      sizeof(VectorScalar) * (dBlock.x*dBlock.y + 2*MAX_COL*dBlock.y) +
      sizeof(size_type) * (MAX_COL*dBlock.y);

#if 0

    std::cout << "Multiply< BlockCrsMatrix< LinearSparse3Tensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << "," << dBlock.z << ")"<< std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  num_spatial_rows(" << num_spatial_rows << ")" << std::endl
              << "  num_stoch_rows(" << num_stoch_rows << ")" << std::endl;
#endif

    //cudaProfilerStart();
    if (A.block.symmetric())
      Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
        ( ApplyKernelSymmetric<MAX_COL>( A, x, y ) );
    else
      Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
        ( ApplyKernelAsymmetric<MAX_COL>( A, x, y ) );
    //cudaProfilerStop();
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_LINEAR_SPARSE_3_TENSOR_HPP */
