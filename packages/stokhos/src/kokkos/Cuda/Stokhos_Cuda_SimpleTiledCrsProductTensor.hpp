// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_CUDA_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP

#include <iostream>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_SimpleTiledCrsProductTensor.hpp"

#include "Stokhos_Cuda_WarpShuffle.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

#if 0

// This version is slow because it requires more shared memory and has
// much more frequent reductions
template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< SimpleTiledCrsProductTensor< TensorScalar, Kokkos::Cuda >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda  execution_space ;
  typedef execution_space::size_type size_type ;

  typedef SimpleTiledCrsProductTensor< TensorScalar, execution_space > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type ;
  typedef Kokkos::View< VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda >  vector_type ;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type m_block_size ;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y,
                       const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), m_block_size(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type max_tile_size = m_A.block.max_jk_tile_size();

      volatile VectorScalar * const sh_x_k =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile VectorScalar * const sh_x_j = sh_x_k+m_block_size*max_tile_size;
      volatile VectorScalar * const sh_A_k = sh_x_j+m_block_size*max_tile_size;
      volatile VectorScalar * const sh_A_j = sh_A_k+m_block_size*max_tile_size;
      volatile VectorScalar * const sh_y   = sh_A_j+m_block_size*max_tile_size;

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / m_block_size;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % m_block_size;
      if (remBlock > 0) ++numBlock;

      // Loop over i tiles
      const size_type n_i_tile = m_A.block.num_i_tiles();
      for (size_type i_tile = 0; i_tile<n_i_tile; ++i_tile) {
        const size_type i_begin = m_A.block.i_begin(i_tile);
        const size_type i_size  = m_A.block.i_size(i_tile);

        // Zero y
        for (size_type i=tid; i<i_size; i+=nid) {
          sh_y[i] = 0.0;
        }

        // Loop over finite element column blocks.
        size_type iBlockEntry = iBlockEntryBeg;
        for (size_type block=0; block<numBlock; ++block,
               iBlockEntry+=m_block_size) {
          const size_type block_size =
            (block == numBlock-1 && remBlock > 0) ? remBlock : m_block_size;

          // Loop over j tiles
          const size_type n_j_tile = m_A.block.num_j_tiles(i_tile);
          for (size_type j_tile = 0; j_tile<n_j_tile; ++j_tile) {
            const size_type j_begin = m_A.block.j_begin(i_tile, j_tile);
            const size_type j_size  = m_A.block.j_size(i_tile, j_tile);

            // Wait for X and A to be used in the previous iteration
            // before reading new values.
            __syncthreads();

            // Coalesced read j-blocks of X and A into shared memory
            for (size_type col=0; col<block_size; ++col) {
              const size_type iBlockColumn =
                m_A.graph.entries(iBlockEntry + col);
              const VectorScalar * const x_j =
                &m_x(j_begin, iBlockColumn);
              const MatrixScalar * const A_j =
                &m_A.values(j_begin, iBlockEntry + col);
              for (size_type j=tid; j<j_size; j+=nid) {
                sh_x_j[col+j*m_block_size] = x_j[j];
                sh_A_j[col+j*m_block_size] = A_j[j];
              }
            }

            // Loop over k tiles
            const size_type n_k_tile = m_A.block.num_k_tiles(i_tile, j_tile);
            for (size_type k_tile = 0; k_tile<n_k_tile; ++k_tile) {
              const size_type k_begin =
                m_A.block.k_begin(i_tile, j_tile, k_tile);
              const size_type k_size  =
                m_A.block.k_size(i_tile, j_tile, k_tile);

              // Wait for X and A to be used in the previous iteration
              // before reading new values.
              __syncthreads();

              // Coalesced read j-blocks of X and A into shared memory
              for (size_type col=0; col<block_size; ++col) {
                const size_type iBlockColumn =
                  m_A.graph.entries(iBlockEntry + col);
                const VectorScalar * const x_k =
                  &m_x(k_begin, iBlockColumn);
                const MatrixScalar * const A_k =
                  &m_A.values(k_begin, iBlockEntry + col);
                for (size_type k=tid; k<k_size; k+=nid) {
                  sh_x_k[col+k*m_block_size] = x_k[k];
                  sh_A_k[col+k*m_block_size] = A_k[k];
                }
              }

              __syncthreads(); // wait for X and A to be read

              // Loop over stochastic rows in this tile
              for (size_type i=threadIdx.y; i<i_size; i+=blockDim.y) {
                VectorScalar s = 0;

                // Product tensor entries which this warp will iterate:
                const size_type lBeg =
                  m_A.block.entry_begin(i_tile, j_tile, k_tile, i);
                const size_type lEnd =
                  m_A.block.entry_end(i_tile, j_tile, k_tile, i);

                // Loop through sparse tensor contributions with
                // coalesced reads.
                for (size_type l=lBeg+threadIdx.x; l<lEnd; l+=blockDim.x) {
                  const size_type kj   = m_A.block.coord( l );
                  const TensorScalar v = m_A.block.value( l );
                  const size_type j    = ( kj & 0x0ffff ) * m_block_size ;
                  const size_type k    = ( kj >> 16     ) * m_block_size ;

                  for ( size_type col = 0; col < block_size; ++col ) {
                    s += v * ( sh_A_j[col+j] * sh_x_k[col+k] +
                               sh_A_k[col+k] * sh_x_j[col+j] );
                  }
                }

                // Reduction of 'y' within 'blockDim.x'
                if (blockDim.x >= 2) s += shfl_down(s, 1, blockDim.x);
                if (blockDim.x >= 4) s += shfl_down(s, 2, blockDim.x);
                if (blockDim.x >= 8) s += shfl_down(s, 4, blockDim.x);
                if (blockDim.x >= 16) s += shfl_down(s, 8, blockDim.x);
                if (blockDim.x >= 32) s += shfl_down(s, 16, blockDim.x);
                if ( threadIdx.x == 0 ) sh_y[i] += s;

              } // i-loop

            } // k-tile loop

          } // j-tile loop

        } // block column loop

        // Wait for all threads to complete the i-tile
        __syncthreads();

        // Store sum for this tile back in global memory
        for (size_type i=tid; i<i_size; i+=nid) {
          m_y( i+i_begin , blockIdx.x ) = sh_y[i];
        }

      } // i-tile loop

    } // operator()

  }; // ProductTensorLoop

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.extent(0) - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type i_tile_size = A.block.max_i_tile_size();
    const size_type jk_tile_size = A.block.max_jk_tile_size();

#ifdef STOKHOS_DEBUG
    const size_type nWarp = 4; // Use fewer warps in debug mode to prevent
                                // launch failures
#else
    const size_type nWarp = 16;
#endif
    const size_type threads_per_row = 16;
    const size_type rows_per_warp =
      Kokkos::Impl::CudaTraits::WarpSize / threads_per_row;
    const dim3 dBlock( threads_per_row , rows_per_warp*nWarp , 1 );
    const dim3 dGrid( row_count , 1 , 1 );

    const size_type shcap =
      Kokkos::Cuda().impl_internal_space_instance()->m_deviceProp.sharedMemPerBlock / 2;
    size_type bs =
      (shcap / sizeof(VectorScalar) - i_tile_size) / (4*jk_tile_size);
    if (bs % 2 == 0) --bs;
    const size_type block_size_max = 31;
    const size_type block_size = std::min(bs, block_size_max);
    // const int block_size = 9;
    const size_type shmem =
      sizeof(VectorScalar) * (4*jk_tile_size*block_size + i_tile_size);

    /*
    ProductTensorLoop kernel( A , x , y, block_size )
    int res;
    cuFuncGetAttribute(&res, CU_FUNC_ATTRIBUTE_NUM_REGS, &kernel.operator());
    */

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tile_size(" << tile_size << ")" << std::endl
              << "  num_tiles(" << num_tiles << ")" << std::endl
             ;
#endif
    //cudaProfilerStart();
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
      ( ProductTensorLoop( A , x , y, block_size ) );
    //cudaProfilerStop();
  }
};

#else

// This is currently the fastest version, but doesn't have fully coalesced
// reads of the sparse tensor, nor coalesced writes of y
template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< SimpleTiledCrsProductTensor< TensorScalar, Kokkos::Cuda >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda  execution_space ;
  typedef execution_space::size_type size_type ;

  typedef SimpleTiledCrsProductTensor< TensorScalar, execution_space > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type ;
  typedef Kokkos::View< VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda >  vector_type ;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type m_block_size ;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y,
                       const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), m_block_size(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type max_tile_size = m_A.block.max_jk_tile_size();

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;
      // const size_type nid2 = nid / 2;
      // const bool lower = tid < nid2;
      // const size_type tid2 = lower ? tid : tid - nid2;

      volatile VectorScalar * const sh_x_k =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile VectorScalar * const sh_x_j = sh_x_k+m_block_size*max_tile_size;
      volatile VectorScalar * const sh_A_k = sh_x_j+m_block_size*max_tile_size;
      volatile VectorScalar * const sh_A_j = sh_A_k+m_block_size*max_tile_size;
      //volatile size_type * const sh_c = (volatile size_type*) (sh_A_j + m_block_size*max_tile_size);
      __shared__ volatile size_type sh_c[32];

      // blockIdx.y == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.y ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.y + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / m_block_size;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % m_block_size;
      if (remBlock > 0) ++numBlock;

      // blockIdx.y == i_tile
      const size_type i_tile = blockIdx.x;
      const size_type i_begin = m_A.block.i_begin(i_tile);
      const size_type i_size  = m_A.block.i_size(i_tile);
      const size_type i1 = threadIdx.y;
      //const size_type i2 = threadIdx.y + blockDim.y;
      if (i1 >= i_size) return;

      VectorScalar s1 = 0;
      //VectorScalar s2 = 0;

      // Loop over finite element column blocks.
      size_type iBlockEntry = iBlockEntryBeg;
      for (size_type block=0; block<numBlock; ++block,
             iBlockEntry+=m_block_size) {
        const size_type block_size =
          (block == numBlock-1 && remBlock > 0) ? remBlock : m_block_size;

        __syncthreads();
        // for (size_type col=tid; col<block_size; col+=nid)
        //   sh_c[col] = m_A.graph.entries(iBlockEntry + col);
        if (tid == 0) {
          for (size_type col=0; col<block_size; ++col)
            sh_c[col] = m_A.graph.entries(iBlockEntry + col);
        }

        // Loop over j tiles
        const size_type n_j_tile = m_A.block.num_j_tiles(i_tile);
        for (size_type j_tile = 0; j_tile<n_j_tile; ++j_tile) {
          const size_type j_begin = m_A.block.j_begin(i_tile, j_tile);
          const size_type j_size  = m_A.block.j_size(i_tile, j_tile);

          // Wait for X and A to be used in the previous iteration
          // before reading new values.
          __syncthreads();

          // Coalesced read j-blocks of X and A into shared memory
          for (size_type col=0; col<block_size; ++col) {
            const VectorScalar * const x_j = &m_x(j_begin, sh_c[col]);
            const MatrixScalar * const A_j = &m_A.values(j_begin, iBlockEntry + col);
            for (size_type j=tid; j<j_size; j+=nid) {
              sh_x_j[col+j*m_block_size] = x_j[j];
              sh_A_j[col+j*m_block_size] = A_j[j];
            }
            // if (lower) {
            //   for (size_type j=tid2; j<j_size; j+=nid2)
            //     sh_x_j[col+j*m_block_size] = x_j[j];
            // }
            // else {
            //   for (size_type j=tid2; j<j_size; j+=nid2)
            //     sh_A_j[col+j*m_block_size] = A_j[j];
            // }
          }

          // Loop over k tiles
          const size_type n_k_tile = m_A.block.num_k_tiles(i_tile, j_tile);
          for (size_type k_tile = 0; k_tile<n_k_tile; ++k_tile) {
            const size_type k_begin =
              m_A.block.k_begin(i_tile, j_tile, k_tile);
            const size_type k_size  =
              m_A.block.k_size(i_tile, j_tile, k_tile);

            // Wait for X and A to be used in the previous iteration
            // before reading new values.
            __syncthreads();

            // Coalesced read j-blocks of X and A into shared memory
            for (size_type col=0; col<block_size; ++col) {
              const VectorScalar * const x_k =
                &m_x(k_begin, sh_c[col]);
              const MatrixScalar * const A_k =
                &m_A.values(k_begin, iBlockEntry + col);
              for (size_type k=tid; k<k_size; k+=nid) {
                sh_x_k[col+k*m_block_size] = x_k[k];
                sh_A_k[col+k*m_block_size] = A_k[k];
              }
              // if (lower) {
              //   for (size_type k=tid2; k<k_size; k+=nid2)
              //     sh_x_k[col+k*m_block_size] = x_k[k];
              // }
              // else {
              //   for (size_type k=tid2; k<k_size; k+=nid2)
              //     sh_A_k[col+k*m_block_size] = A_k[k];
              // }
            }

            __syncthreads(); // wait for X and A to be read

            // Product tensor entries which this warp will iterate:
            size_type lBeg =
              m_A.block.entry_begin(i_tile, j_tile, k_tile, i1);
            size_type lEnd =
              m_A.block.entry_end(i_tile, j_tile, k_tile, i1);

            // Loop through sparse tensor contributions with
            // coalesced reads (these aren't fully coalesced unless
            // blockDim.x >= 16 for double precision.  We need to reorder
            // the tensor entries for fewer threads per row).
            for (size_type l=lBeg+threadIdx.x; l<lEnd; l+=blockDim.x) {
              const size_type kj   = m_A.block.coord( l );
              const TensorScalar v = m_A.block.value( l );
              const size_type j    = ( kj & 0x0ffff ) * m_block_size ;
              const size_type k    = ( kj >> 16     ) * m_block_size ;

              for ( size_type col = 0; col < block_size; ++col ) {
                s1 += v * ( sh_A_j[col+j] * sh_x_k[col+k] +
                            sh_A_k[col+k] * sh_x_j[col+j] );
              }
            }

            // if (i2 < i_size) {
            //   // Product tensor entries which this warp will iterate:
            //   size_type lBeg =
            //     m_A.block.entry_begin(i_tile, j_tile, k_tile, i2);
            //   size_type lEnd =
            //     m_A.block.entry_end(i_tile, j_tile, k_tile, i2);

            //   // Loop through sparse tensor contributions with
            //   // coalesced reads.
            //   for (size_type l=lBeg+threadIdx.x; l<lEnd; l+=blockDim.x) {
            //     const size_type kj   = m_A.block.coord( l );
            //     const TensorScalar v = m_A.block.value( l );
            //     const size_type j    = ( kj & 0x0ffff ) * m_block_size ;
            //     const size_type k    = ( kj >> 16     ) * m_block_size ;

            //     for ( size_type col = 0; col < block_size; ++col ) {
            //       s2 += v * ( sh_A_j[col+j] * sh_x_k[col+k] +
            //                   sh_A_k[col+k] * sh_x_j[col+j] );
            //     }
            //   }
            // }

          } // k-tile loop

        } // j-tile loop

      } // block column loop

      // Wait for all threads to complete the i-tile
      //__syncthreads();

      // Reduction of 'y' within 'blockDim.x'
      if (blockDim.x >= 2) s1 += shfl_down(s1, 1, blockDim.x);
      if (blockDim.x >= 4) s1 += shfl_down(s1, 2, blockDim.x);
      if (blockDim.x >= 8) s1 += shfl_down(s1, 4, blockDim.x);
      if (blockDim.x >= 16) s1 += shfl_down(s1, 8, blockDim.x);
      if (blockDim.x >= 32) s1 += shfl_down(s1, 16, blockDim.x);
      if ( threadIdx.x == 0 ) m_y( i1+i_begin , blockIdx.y ) = s1;

      // if (i2 < i_size) {
      //   if (blockDim.x >= 2) s2 += shfl_down(s2, 1, blockDim.x);
      //   if (blockDim.x >= 4) s2 += shfl_down(s2, 2, blockDim.x);
      //   if (blockDim.x >= 8) s2 += shfl_down(s2, 4, blockDim.x);
      //   if (blockDim.x >= 16) s2 += shfl_down(s2, 8, blockDim.x);
      //   if (blockDim.x >= 32) s2 += shfl_down(s2, 16, blockDim.x);
      //   if ( threadIdx.x == 0 ) m_y( i2+i_begin , blockIdx.y ) = s2;
      // }

    } // operator()

  }; // ProductTensorLoop

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.extent(0) - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tile_size = A.block.max_jk_tile_size();
    const size_type n_i_tile = A.block.num_i_tiles();

#ifdef STOKHOS_DEBUG
    const size_type nWarp = 2; // Use fewer warps in debug mode to prevent
                                // launch failures
#else
    const size_type nWarp = 16;
#endif
    const size_type threads_per_row = 4;
    const size_type rows_per_warp =
      Kokkos::Impl::CudaTraits::WarpSize / threads_per_row;
    const dim3 dBlock( threads_per_row , rows_per_warp*nWarp , 1 );
    const dim3 dGrid( n_i_tile , row_count , 1 );

    const size_type shcap =
      Kokkos::Cuda().impl_internal_space_instance()->m_deviceProp.sharedMemPerBlock / 2;
    size_type bs = ((shcap / sizeof(VectorScalar)) / tile_size) / 4;
    if (bs % 2 == 0) --bs;
    const size_type block_size_max = 31;
    const size_type block_size = std::min(bs, block_size_max);
    // const int block_size = 9;
    const size_type shmem = sizeof(VectorScalar)*(4*block_size*tile_size); // + sizeof(size_type)*block_size_max;

    // Try sorting rows based on number of non-zeros, split threads between
    // A, x reads, do the diagonals properly
    // Splitting threads between A & x reads doesn't seem to improve things

    /*
    ProductTensorLoop kernel( A , x , y, block_size )
    int res;
    cuFuncGetAttribute(&res, CU_FUNC_ATTRIBUTE_NUM_REGS, &kernel.operator());
    */

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tile_size(" << tile_size << ")" << std::endl
              << "  num_i_tiles(" << n_i_tile << ")" << std::endl
             ;
#endif
    //cudaProfilerStart();
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
      ( ProductTensorLoop( A , x , y, block_size ) );
    //cudaProfilerStop();
  }
};

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_TILED_CRS_PRODUCT_TENSOR_HPP */
