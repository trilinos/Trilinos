// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_CUDA_CRS_PRODUCT_TENSOR_HPP

#include <iostream>

#include "Kokkos_Cuda.hpp"
#include "Cuda/Kokkos_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_CrsProductTensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

#if 1

//----------------------------------------------------------------------------

template< typename TensorScalar,
          typename MatrixScalar,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< TensorScalar, Kokkos::Cuda >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  DefaultSparseMatOps >
{
public:

  typedef Kokkos::Cuda                    device_type;
  typedef device_type::size_type  size_type;

  typedef CrsProductTensor< TensorScalar, device_type >       tensor_type;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type;
  typedef Kokkos::View< VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda >           vector_type;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type BlockSize;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y,
                       const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), BlockSize(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      VectorScalar * const sh_x =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_A = sh_x + BlockSize*dim;
      VectorScalar * const sh_y = sh_A + BlockSize*dim;
      volatile VectorScalar * const sh_t = sh_y + dim;

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / BlockSize;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % BlockSize;
      if (remBlock > 0) ++numBlock;

      // Zero y
      for ( size_type i = tid; i < dim; i += nid ) {
        sh_y[i] = 0.0;
      }

      // Loop over columns in the discrete (finite element) system.
      size_type iBlockEntry = iBlockEntryBeg;
      for ( size_type block = 0; block < numBlock;
            ++block, iBlockEntry += BlockSize) {
        const size_type block_size =
          (block == numBlock-1 && remBlock > 0) ? remBlock : BlockSize;

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
              y += v * ( sh_A[col+j] * sh_x[col+k] + sh_A[col+k] * sh_x[col+j] );
            }

          }

          // Reduction of 'y' within 'blockDim.x'
          sh_t[ tid ] = y;
          if (threadIdx.x+16 < blockDim.x) sh_t[tid] += sh_t[tid+16];
          if (threadIdx.x+ 8 < blockDim.x) sh_t[tid] += sh_t[tid+ 8];
          if (threadIdx.x+ 4 < blockDim.x) sh_t[tid] += sh_t[tid+ 4];
          if (threadIdx.x+ 2 < blockDim.x) sh_t[tid] += sh_t[tid+ 2];
          if (threadIdx.x+ 1 < blockDim.x) sh_t[tid] += sh_t[tid+ 1];
          if (threadIdx.x == 0) sh_y[i] += sh_t[tid];

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

  //------------------------------------

  static void apply( const matrix_type & A,
                     const vector_type & x,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.dimension_0() - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tensor_align = tensor_dimension;

#ifdef STOKHOS_DEBUG
    const size_type nWarp = 16; // Use fewer warps in debug mode to prevent
                                // launch failures
#else
     const size_type nWarp = 20;
#endif
    const dim3 dBlock( Kokkos::Impl::CudaTraits::WarpSize, nWarp, 1 );
    const dim3 dGrid( row_count, 1, 1 );

    // Use at most half of shared memory to get 2 blocks per SMP
    const size_type shcap =
      Kokkos::Impl::CudaTraits::SharedMemoryCapacity / 2;
    int block_size = ((shcap / sizeof(VectorScalar) - dBlock.x*dBlock.y) / tensor_align - 1) / 2;
    block_size = std::min( block_size, 9 );
    if (block_size % 2 == 0)
      --block_size;
    // const int block_size = 9;
    const size_type shmem =
      sizeof(VectorScalar) * ((2*block_size+1) * tensor_align + dBlock.x*dBlock.y);

#if 0

    const size_type mega = 1024*1024;
    const size_type giga = 1024*mega;
    const size_type fem_nnz = A.values.dimension_1();
    const size_type tensor_nnz = A.block.num_non_zeros();
    const size_type A_reads =
      ((static_cast<double>(fem_nnz) * tensor_dimension)*sizeof(MatrixScalar))/mega;
    const size_type x_reads =
      ((static_cast<double>(fem_nnz) * tensor_dimension)*sizeof(VectorScalar))/mega;
    const size_type y_writes =
      ((static_cast<double>(row_count) * tensor_dimension)*sizeof(VectorScalar))/mega;
    const size_type tensor_reads =
      ((static_cast<double>(fem_nnz) / block_size)*tensor_nnz*(sizeof(TensorScalar)+2*sizeof(size_type)))/giga;

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  fem_nnz(" << fem_nnz << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_entry_max(" << A.block.entry_maximum() << ")" << std::endl
              << "  A+x+y reads/writes(" << A_reads+x_reads+y_writes << " MB)" << std::endl
              << "  tensor reads(" << tensor_reads << " GB)" << std::endl
             ;
#endif
    //cudaProfilerStart();
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
      ( ProductTensorLoop( A, x, y, block_size ) );
    //cudaProfilerStop();
  }
};

#else

//----------------------------------------------------------------------------

/* Working on splitting warps between stochastic rows.  Requires reading
 * tensor values into shared memory, which is getting complicated.  Doesn't
 * work yet.
 */

template< typename TensorScalar,
          typename MatrixScalar,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< TensorScalar, Kokkos::Cuda >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  DefaultSparseMatOps >
{
public:

  typedef Kokkos::Cuda                    device_type;
  typedef device_type::size_type  size_type;

  typedef CrsProductTensor< TensorScalar, device_type >       tensor_type;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type;
  typedef Kokkos::View< VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda >           vector_type;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type BlockSize;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y,
                       const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), BlockSize(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      const size_type nid = blockDim.x * blockDim.y * blockDim.z;
      const size_type tid =
        threadIdx.x + blockDim.x * ( threadIdx.y + blockDim.y * threadIdx.z );
      const size_type row_lane = threadIdx.y + blockDim.y * threadIdx.z;
      const size_type nrow = blockDim.y * blockDim.z;
      const size_type warp_size = blockDim.x * blockDim.y;
      const size_type thread_lane = threadIdx.x + blockDim.x * threadIdx.y;

      VectorScalar * const sh_x =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_A = sh_x + BlockSize*dim;
      VectorScalar * const sh_y = sh_A + BlockSize*dim;
      volatile VectorScalar * const sh_t = sh_y + dim;
      size_type * const sh_j = (size_type*) (sh_t + nid);
      size_type * const sh_k = (size_type*) (sh_j + nid);

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / BlockSize;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % BlockSize;
      if (remBlock > 0) ++numBlock;

      // Zero y
      for ( size_type i = tid; i < dim; i += nid ) {
        sh_y[i] = 0.0;
      }

      // Loop over columns in the discrete (finite element) system.
      size_type iBlockEntry = iBlockEntryBeg;
      for ( size_type block = 0; block < numBlock;
            ++block, iBlockEntry += BlockSize) {
        const size_type block_size =
          (block == numBlock-1 && remBlock > 0) ? remBlock : BlockSize;

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
        for (size_type bi = blockDim.y*threadIdx.z; bi < dim; bi += nrow) {

          // Read tensor entries
          for (size_type irow = 0; irow < blockDim.y; ++irow) {
            const size_type i = irow + bi;
            if (i < dim) {

              /*
               * This needs work!  We don't want to read the whole row
               * into shared memory.
               */

              // Product tensor entries which this warp will iterate:
              const size_type lBeg = m_A.block.entry_begin( bi+irow );
              const size_type lEnd = m_A.block.entry_end(   bi+irow );

              // Loop through sparse tensor contributions with coalesced reads.
              for (size_type l = lBeg+thread_lane; l < lEnd; l += warp_size) {

                // Read 'blockDim.x' entries from the tensor
                sh_j[thread_lane+] = m_A.block.coord(l,0);
                sh_k[tid] = m_A.block.coord(l,1);
                sh_t[tid] = m_A.block.value(l);

              }

            }

          }
          VectorScalar y = 0;
          const size_type i = bi + threadIdx.y;
          if (i < dim) {

            // Product tensor entries which this warp will iterate:
            const size_type lBeg = m_A.block.entry_begin( i );
            const size_type lEnd = m_A.block.entry_end(   i );

            // Loop through sparse tensor contributions with coalesced reads.
            for (size_type l = lBeg; l < lEnd; l += warp_size) {

              const int mm = (l+warp_size<lEnd) ? warp_size : lEnd-l;
              for (size_type m=threadIdx.x; m<mm; m+=blockDim.x) {

                // const size_type j = m_A.block.coord(l+m,0);
                // const size_type k = m_A.block.coord(l+m,1);
                // const TensorScalar v = m_A.block.value(l+m);
                const size_type j = sh_j[m+warp_size*threadIdx.z];
                const size_type k = sh_k[m+warp_size*threadIdx.z];
                const TensorScalar v = sh_t[m+warp_size*threadIdx.z];

                for (size_type col = 0; col < block_size; ++col) {
                  const size_type jj = col + j * BlockSize;
                  const size_type kk = col + k * BlockSize;
                  y += v * ( sh_A[jj] * sh_x[kk] + sh_A[kk] * sh_x[jj] );
                }

              }

            }

          }

          /*
          for (size_type l = lBeg + threadIdx.x; l < lEnd; l += blockDim.x) {

            // Read 'blockDim.x' entries from the tensor
            const size_type j = m_A.block.coord( l, 0 ); // coalesced read
            const size_type k = m_A.block.coord( l, 1 ); // coalesced read
            const TensorScalar v = m_A.block.value(l);   // coalesced read

            for (size_type col = 0; col < block_size; ++col) {
              const size_type jj = col + j * BlockSize;
              const size_type kk = col + k * BlockSize;
              y += v * ( sh_A[jj] * sh_x[kk] + sh_A[kk] * sh_x[jj] );
            }

          }
          */

          // Reduction of 'y' within 'blockDim.x'
          sh_t[ tid ] = y;
          if (threadIdx.x + 16 < blockDim.x) sh_t[tid] += sh_t[tid+16];
          if (threadIdx.x +  8 < blockDim.x) sh_t[tid] += sh_t[tid+ 8];
          if (threadIdx.x +  4 < blockDim.x) sh_t[tid] += sh_t[tid+ 4];
          if (threadIdx.x +  2 < blockDim.x) sh_t[tid] += sh_t[tid+ 2];
          if (threadIdx.x +  1 < blockDim.x) sh_t[tid] += sh_t[tid+ 1];
          if (threadIdx.x == 0) sh_y[i] += sh_t[tid];

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

  //------------------------------------

  static void apply( const matrix_type & A,
                     const vector_type & x,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.dimension_0() - 1;
    const size_type tensor_dimension = A.block.dimension();
    // const size_type rem = tensor_dimension % 32;
    // const size_type tensor_align =
    //   rem > 0 ? tensor_dimension + rem : tensor_dimension;
    const size_type tensor_align = tensor_dimension;

    const size_type warp_size = Kokkos::Impl::CudaTraits::WarpSize;
    const size_type row_size = 4;
    const size_type nWarp = 16;
    const dim3 dBlock( warp_size / row_size, row_size, nWarp );
    const dim3 dGrid( row_count, 1, 1 );

    // Use at most half of shared memory to get 2 blocks per SMP
    const size_type shcap =
      Kokkos::Impl::CudaTraits::SharedMemoryCapacity / 2;
    int block_size = (((shcap - sizeof(size_type) * 2 * dBlock.x*dBlock.y*dBlock.z) / sizeof(VectorScalar) - dBlock.x*dBlock.y*dBlock.z) / tensor_align - 1) / 2;
    block_size = std::min( block_size, 9 );
    if (block_size % 2 == 0)
      --block_size;
    //const int block_size = 9;
    const size_type shmem =
      sizeof(VectorScalar) * ((2*block_size+1) * tensor_align + dBlock.x*dBlock.y*dBlock.z) +
      sizeof(size_type) * 2 * dBlock.x*dBlock.y*dBlock.z;

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << "," << dBlock.z
              << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_entry_max(" << A.block.entry_maximum() << ")" << std::endl
              ;
#endif
    //cudaProfilerStart();
    Kokkos:: Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
      ( ProductTensorLoop( A, x, y, block_size ) );
    //cudaProfilerStop();
  }
};

#endif
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_CRS_PRODUCT_TENSOR_HPP */
