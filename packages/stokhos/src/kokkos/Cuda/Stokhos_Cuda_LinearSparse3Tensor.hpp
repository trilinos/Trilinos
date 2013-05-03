// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_FLAT_SPARSE_3_TENSOR_HPP
#define STOKHOS_CUDA_FLAT_SPARSE_3_TENSOR_HPP

#include <iostream>

#include "KokkosArray_Cuda.hpp"
#include "Cuda/KokkosArray_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_LinearSparse3Tensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar ,
          int BlockSize >
class Multiply<
  BlockCrsMatrix< LinearSparse3Tensor< TensorScalar, KokkosArray::Cuda , BlockSize >,
                  MatrixScalar, KokkosArray::Cuda >,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Cuda device_type ;
  typedef device_type::size_type size_type ;

  typedef LinearSparse3Tensor< TensorScalar , device_type > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef KokkosArray::View< VectorScalar** ,
                             KokkosArray::LayoutLeft ,
                             KokkosArray::Cuda > vector_type ;

  class ApplyKernel {
  public:

    const matrix_type m_A ;
    const vector_type m_x ;
    const vector_type m_y ;

    ApplyKernel( const matrix_type & A ,
                 const vector_type & x ,
                 const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      const size_type nid = blockDim.x * blockDim.y ;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y ;

      VectorScalar * const sh_t =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      volatile VectorScalar * const sh_y0 = sh_t + tid;
      VectorScalar * const sh_y = sh_t + nid;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];

      // Read tensor entries
      const TensorScalar c0 = m_A.block.value(0);
      const TensorScalar c1 = m_A.block.value(1);

      // Zero y
      for ( size_type i = tid ; i < dim ; i += nid ) {
        sh_y[i] = 0.0;
      }
      sh_y0[0] = 0.0;
      VectorScalar y0 = 0.0;

      // Loop over columns in the discrete (finite element) system.
      for ( size_type iBlockEntry = iBlockEntryBeg ;
            iBlockEntry != iBlockEntryEnd ; ++iBlockEntry ) {

        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry  );

        const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

        const VectorScalar a0 = A[0];
        const MatrixScalar x0 = x[0];
        y0 -= 2.0*c0*a0*x0;

        // Loop over rows of stochastic block
        for ( size_type i = tid ; i < dim ; i += nid) {
          const MatrixScalar ai = A[i];
          const VectorScalar xi = x[i];
          sh_y[i] += c1*(a0*xi + ai*x0);
          y0  += c1*ai*xi;
        }

      }

      sh_y0[ tid ] = y0 ;
      __syncthreads();

      // Reduction of 'y0' across warps
      size_type warp = blockDim.y >> 1;
      while (warp >= 1) {
        if ( threadIdx.y + warp < blockDim.y ) 
          sh_y0[tid] += sh_y0[tid+warp*blockDim.x];
        __syncthreads();
        warp >>= 1;
      }

      // Reduction of 'y0' within warp
      if ( threadIdx.x + 16 < blockDim.x ) sh_y0[tid] += sh_y0[tid+16];
      if ( threadIdx.x +  8 < blockDim.x ) sh_y0[tid] += sh_y0[tid+ 8];
      if ( threadIdx.x +  4 < blockDim.x ) sh_y0[tid] += sh_y0[tid+ 4];
      if ( threadIdx.x +  2 < blockDim.x ) sh_y0[tid] += sh_y0[tid+ 2];
      if ( threadIdx.x +  1 < blockDim.x ) sh_y0[tid] += sh_y0[tid+ 1];

      if ( threadIdx.x == 0 )
        sh_y[0] += sh_y0[0];

      __syncthreads();

      // Store result back in global memory
      for ( size_type i = tid ; i < dim ; i += nid ) {
        m_y( i , blockIdx.x ) = sh_y[ i ] ;
      }
    }
  };


  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type num_spatial_rows = A.graph.row_map.dimension_0() - 1 ;
    const size_type num_stoch_rows = A.block.dimension();

    const size_type warp_size = KokkosArray::Impl::CudaTraits::WarpSize;
    // const size_type max_warp =
    //   KokkosArray::Impl::cuda_internal_maximum_warp_count();
    // const size_type nWarp = std::min(num_stoch_rows / warp_size, max_warp);
    const size_type nWarp = 16;
    const dim3 dBlock( warp_size , nWarp , 1 );
    const dim3 dGrid( num_spatial_rows , 1 , 1 );

    const size_type shmem =
      sizeof(VectorScalar) * (num_stoch_rows + dBlock.x*dBlock.y);

#if 1

    std::cout << "Multiply< BlockCrsMatrix< LinearSparse3Tensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  num_spatial_rows(" << num_spatial_rows << ")" << std::endl
              << "  num_stoch_rows(" << num_stoch_rows << ")" << std::endl;
#endif

    //cudaProfilerStart();
    KokkosArray:: Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
      ( ApplyKernel( A , x , y ) );
    //cudaProfilerStop();
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_LINEAR_SPARSE_3_TENSOR_HPP */
