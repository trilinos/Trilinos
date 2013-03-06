/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CUDA_CRSPRODUCTTENSORLEGENDRE_HPP
#define KOKKOSARRAY_CUDA_CRSPRODUCTTENSORLEGENDRE_HPP

#include <stdio.h>

#include <utility>
#include <sstream>
#include <stdexcept>
#include <KokkosArray_BlockCrsMatrix.hpp>
#include <KokkosArray_CrsProductTensorLegendre.hpp>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_NONE                0  /* verified */
#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_SHARED_A_X          1  /* verified */
#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_SHARED_MULTI_A_X    2  /* verified */
#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_SHARED_BLOCK_A_X    3  /* verified */
#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_GLOBAL_BLOCK_A_X    4  /* verified */
#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_VARIABLE_BLOCK_A_X  5  /* verified */

#define MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT  2

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT == MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_NONE

namespace KokkosArray {
namespace Impl {

/* Non-specialized algorithm to compare for correctness */

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
private:

  typedef BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                          MatrixScalar , Cuda > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda > vector_type ;

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;

public:

  typedef Cuda                             device_type ;
  typedef typename device_type::size_type  size_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type iy ) const
  {
    // Compute spatial row 'iy'

    const size_type tensor_dim     = m_A.block.dimension();
    const size_type iBlockEntryEnd = m_A.graph.row_map[ iy + 1 ];
          size_type iBlockEntry    = m_A.graph.row_map[ iy ];

    VectorScalar * const y = & m_y( 0 , iy );

    for ( size_type iyInner = 0 ; iyInner < tensor_dim ; ++iyInner ) {
      y[iyInner] = 0 ;
    }

    for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {
      const VectorScalar * const x = & m_x( 0 , m_A.graph.entries( iBlockEntry ) );
      const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

      m_A.block.multiply_add( A , x , y );
    }
  }

  Multiply( const matrix_type & arg_A ,
            const vector_type & arg_x ,
            const vector_type & arg_y )
  : m_A( arg_A )
  , m_x( arg_x )
  , m_y( arg_y )
  {}

  void run() const
  {
    parallel_for( m_y.dimension_1() , *this );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT == MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_SHARED_MULTI_A_X

namespace KokkosArray {
namespace Impl {

/** \brief
 *
 *  finite_elem_row = blockIdx.x ;
 *
 *  Read all of stochastic block coefficients for A and x into shared memory
 *  and perform random accesses into those arrays.
 *
 *  Re-read (coalesced) tensor product for each finite element column.
 */

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensorLegendre< MatrixScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type >  matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >                vector_type ;

private:

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;
  mutable size_type m_block_size ;

public:

  __device__
  void operator()(void) const
  {
    // Number of bases in the stochastic system and padding to warp size.
    const size_type tensor_dim = m_A.block.dimension();
    const size_type tensor_dim_align = tensor_dim ;

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_y[ tensor_dim_align ]
    //   sh_x[ tensor_dim_align * m_block_size ]
    //   sh_A[ tensor_dim_align * m_block_size ]
    //   sh_offset[ 2 * tensor_dim_align + 1 ]

    volatile VectorScalar * const sh_work = kokkos_impl_cuda_shared_memory<VectorScalar>();
    VectorScalar * const sh_y = (VectorScalar *)( sh_work + blockDim.x * blockDim.y );
    VectorScalar * const sh_x = (VectorScalar *)( sh_y + tensor_dim_align );
    MatrixScalar * const sh_A = (MatrixScalar *)( sh_x + tensor_dim_align * m_block_size );
    //unsigned     * const sh_offset = (unsigned *)( sh_A + tensor_dim_align * m_block_size );

    // Block size and thread id within the entire block:

    const size_type nid = CudaTraits::WarpSize * blockDim.y ;
    const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

    // { // Load the tensor offsets into shared memory
    //   const size_type n = 2 * tensor_dim + 1 ;
    //   for ( size_type i = tid ; i < n ; ++i ) {
    //     sh_offset[i] = m_A.block.m_entry_offset(i);
    //   }
    // }

    // blockIdx.x == row in the deterministic (finite element) system
    const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
    const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];

    size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / m_block_size;
    const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % m_block_size;
    if (remBlock > 0) ++numBlock;

    // Zero y
    for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
      sh_y[i] = 0.0;
    }

     // Loop over columns in the discrete (finite element) system.
    size_type iBlockEntry = iBlockEntryBeg ;
    for ( size_type block = 0 ; block < numBlock ; ++block ) {
      const size_type block_size = 
	(block == numBlock-1 && remBlock > 0) ? remBlock : m_block_size;

      // Coalesced read blocks of X and A into shared memory
      for ( size_type iBlock = 0; iBlock < block_size ; ++iBlock ) {

	const size_type iBlockColumn = 
	  m_A.graph.entries( iBlockEntry + iBlock );

	const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry + iBlock );

        // Wait for X and A to be used in the previous iteration before reading new values.
        __syncthreads();

	// Coalesced read by the whole block from global memory:
        for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
          sh_x[iBlock+i*block_size] = x[i] ; // m_x(        i , iBlockColumn );
          sh_A[iBlock+i*block_size] = A[i] ; // m_A.values( i , iBlockEntry );
        }

      }

      __syncthreads();
      // Wait for X and A to be read before using these values in the next iteration.

      // Each warp computes a row of the stochastic block 
      // for coalesced reads and no need for explicit synchronization.

      for ( size_type iyInner = threadIdx.y ; iyInner < tensor_dim ; iyInner += blockDim.y ) {

        VectorScalar y = 0 ;

        // Product tensor entries which this warp will iterate:
        //
        // const size_type iBeg        = sh_offset[ 2 * iyInner ];
        // const size_type iBegOffDiag = sh_offset[ 2 * iyInner + 1 ];
        // const size_type iEnd        = sh_offset[ 2 * iyInner + 2 ];

	const size_type iBeg        = m_A.block.m_entry_offset( 2 * iyInner );
        const size_type iBegOffDiag = m_A.block.m_entry_offset( 2 * iyInner + 1 );
        const size_type iEnd        = m_A.block.m_entry_offset( 2 * iyInner + 2 );

        // Loop through sparse tensor diagonal contributions:

        for ( size_type i = threadIdx.x + iBeg ; i < iBegOffDiag ; i += blockDim.x ) {
          const unsigned j = m_A.block.m_coordinate(i);
	  const MatrixScalar v = m_A.block.m_value(i);
	  for ( size_type iBlock = 0; iBlock < block_size ; ++iBlock ) {
	    const size_type jj = iBlock+j*block_size;
	    y += v * sh_A[jj] * sh_x[jj] ;
	  }
        }

        // Loop through sparse tensor off-diagonal contributions:

        for ( size_type i = threadIdx.x + iBegOffDiag ; i < iEnd ; i += blockDim.x ) {
          const unsigned kj = m_A.block.m_coordinate(i);
          const unsigned j  = kj & 0x0ffff ;
          const unsigned k  = kj >> 16 ;
	  const MatrixScalar v = m_A.block.m_value(i);
	  for ( size_type iBlock = 0; iBlock < block_size ; ++iBlock ) {
	    const size_type jj = iBlock+j*block_size;
	    const size_type kk = iBlock+k*block_size;
	    y += v * ( sh_A[jj] * sh_x[kk] + sh_A[kk] * sh_x[jj] );
	  }
        }

        // Reduction of 'y' within the warp

        sh_work[ tid ] = y ;

        if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+16];
        if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 8];
        if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 4];
        if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 2];
        if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 1];

        if ( 0 == threadIdx.x ) { // One thread in the warp saves to shared memory
          sh_y[ iyInner ] += sh_work[ tid ];
        }
      }

      iBlockEntry += block_size;
    }

    __syncthreads();

    // Coalesced write by the whole block to global memory
    for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
      m_y( i , blockIdx.x ) = sh_y[i];
    }
  }

  //------------------------------------

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
  : m_A( A ), m_x( x ), m_y( y ) {}

  void run() const
  {
    const size_type row_count        = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dimension = m_A.block.dimension();
    const size_type tensor_dim_align = tensor_dimension ;

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_y[ tensor_dim_align ]
    //   sh_x[ tensor_dim_align * block_size ]
    //   sh_A[ tensor_dim_align * block_size ]
    //   sh_offset[ 2 * tensor_dim_align + 1 ]

    const size_type nWarp = 16 ; 
    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );
    const dim3 dGrid( row_count , 1 , 1 );

    // Use at most half of shared memory to get 2 blocks per SMP
    const size_type shcap = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity / 2;
    //- sizeof(unsigned int) * ( 2 * tensor_dim_align + 1 );
    m_block_size = ((shcap / sizeof(VectorScalar) - dBlock.x*dBlock.y) / tensor_dim_align - 1) / 2;
    m_block_size = std::min( m_block_size, size_type(9) );
    if (m_block_size % 2 == 0)
      --m_block_size;
    //const int m_block_size = 5;
    const size_type shmem = 
      sizeof(VectorScalar) * ((2*m_block_size+1) * tensor_dim_align + dBlock.x*dBlock.y);
      //sizeof(unsigned int) * ( 2 * tensor_dim_align + 1 );

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensorLegendre ... > >"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
	      << "  block_size(" << m_block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_dim_align(" << tensor_dim_align << ")" << std::endl
              ;
#endif

    cuda_parallel_launch_local_memory< Multiply ><<< dGrid , dBlock , shmem >>>( *this );
    // Impl::CudaParallelLaunch< Multiply >( *this, dGrid , dBlock , shmem );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT == MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_SHARED_A_X

namespace KokkosArray {
namespace Impl {

/** \brief
 *
 *  finite_elem_row = blockIdx.x ;
 *
 *  Read all of stochastic block coefficients for A and x into shared memory
 *  and perform random accesses into those arrays.
 *
 *  Re-read (coalesced) tensor product for each finite element column.
 */

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensorLegendre< MatrixScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type >  matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >                vector_type ;

private:

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;

public:

  __device__
  void operator()(void) const
  {
    // Number of bases in the stochastic system and padding to warp size.
    const size_type tensor_dim = m_A.block.dimension();
    const size_type tensor_dim_align = tensor_dim ;

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_y[ tensor_dim_align ]
    //   sh_x[ tensor_dim_align ]
    //   sh_A[ tensor_dim_align ]
    //   sh_offset[ 2 * tensor_dim_align + 1 ]

    volatile VectorScalar * const sh_work = kokkos_impl_cuda_shared_memory<VectorScalar>();
    VectorScalar * const sh_y = (VectorScalar *)( sh_work + blockDim.x * blockDim.y );
    VectorScalar * const sh_x = (VectorScalar *)( sh_y + tensor_dim_align );
    MatrixScalar * const sh_A = (MatrixScalar *)( sh_x + tensor_dim_align );
    unsigned     * const sh_offset = (unsigned *)( sh_A + tensor_dim_align );

    // Block size and thread id within the entire block:

    const size_type nid = CudaTraits::WarpSize * blockDim.y ;
    const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

    { // Zero out values of 'y' for accumulation across discrete columns
      for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
        sh_y[i] = 0 ;
      }
    }

    { // Load the tensor offsets into shared memory
      const size_type n = 2 * tensor_dim + 1 ;
      for ( size_type i = tid ; i < n ; ++i ) {
        sh_offset[i] = m_A.block.m_entry_offset(i);
      }
    }

    // Loop over columns in the discrete (finite element) system.

    // blockIdx.x == row in the deterministic (finite element) system
    const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
          size_type iBlockEntry    = m_A.graph.row_map[ blockIdx.x ];

    for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

      { // Read stochastic coefficients from the vector and matrix into shared memory.
        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );

        const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

        // Wait for X and A to be used in the previous iteration before reading new values.
        __syncthreads();

        // Coalesced read by the whole block from global memory:
        for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
          sh_x[i] = x[i] ; // m_x(        i , iBlockColumn );
          sh_A[i] = A[i] ; // m_A.values( i , iBlockEntry );
        }

        __syncthreads();
        // Wait for X and A to be read before using these values in the next iteration.
      }

      // Each warp computes a row of the stochastic block 
      // for coalesced reads and no need for explicit synchronization.

      for ( size_type iyInner = threadIdx.y ; iyInner < tensor_dim ; iyInner += blockDim.y ) {

        VectorScalar y = 0 ;

        // Product tensor entries which this warp will iterate:
        //
        const size_type iBeg        = sh_offset[ 2 * iyInner ];
        const size_type iBegOffDiag = sh_offset[ 2 * iyInner + 1 ];
        const size_type iEnd        = sh_offset[ 2 * iyInner + 2 ];

        // Loop through sparse tensor diagonal contributions:

        for ( size_type i = iBeg + threadIdx.x ; i < iBegOffDiag ; i += blockDim.x ) {
          const unsigned j = m_A.block.m_coordinate(i);
          y += m_A.block.m_value(i) * sh_A[j] * sh_x[j] ;
        }

        // Loop through sparse tensor off-diagonal contributions:

        for ( size_type i = iBegOffDiag + threadIdx.x ; i < iEnd ; i += blockDim.x ) {
          const unsigned kj = m_A.block.m_coordinate(i);
          const unsigned j  = kj & 0x0ffff ;
          const unsigned k  = kj >> 16 ;
          y += m_A.block.m_value(i) * ( sh_A[j] * sh_x[k] + sh_A[k] * sh_x[j] );
        }

        // Reduction of 'y' within the warp

        sh_work[ tid ] = y ;

        if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+16];
        if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 8];
        if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 4];
        if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 2];
        if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 1];

        if ( 0 == threadIdx.x ) { // One thread in the warp saves to shared memory
          sh_y[ iyInner ] += sh_work[ tid ];
        }
      }
    }

    __syncthreads();

    // Coalesced write by the whole block to global memory
    for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
      m_y( i , blockIdx.x ) = sh_y[i];
    }
  }

  //------------------------------------

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
  : m_A( A ), m_x( x ), m_y( y ) {}

  void run() const
  {
    const size_type row_count        = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dimension = m_A.block.dimension();
    const size_type tensor_dim_align = tensor_dimension ;

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_y[ tensor_dim_align ]
    //   sh_x[ tensor_dim_align ]
    //   sh_A[ tensor_dim_align ]
    //   sh_offset[ 2 * tensor_dim_align + 1 ]

    const int nsh_warp = sizeof(VectorScalar) * CudaTraits::WarpSize ; // sh_work

    const int nsh_base = sizeof(VectorScalar) * tensor_dim_align             // sh_y
                       + sizeof(VectorScalar) * tensor_dim_align             // sh_x
                       + sizeof(MatrixScalar) * tensor_dim_align             // sh_A
                       + sizeof(unsigned int) * ( 2 * tensor_dim_align + 1 ) // sh_offset
                       ;

    enum { ShCapacity = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity };

    int nWarp = ( ShCapacity - nsh_base ) / nsh_warp ;

    if ( nWarp < 1 ) {
      std::ostringstream msg ;
      msg << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply "
          << " FAILED to meet shared memory need: "
          << nsh_base << " + " << nsh_warp
          << " * nWarp <= " << ShCapacity ;
      throw std::runtime_error( msg.str() );
    }

    nWarp = std::min( nWarp , (int) cuda_internal_maximum_warp_count() );
    nWarp = std::min( nWarp , (int) tensor_dimension );

    const size_type shmem = nsh_base + nsh_warp * nWarp ;

    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );

    const dim3 dGrid( row_count , 1 , 1 );

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensorLegendre ... > >"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem << ")" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_dim_align(" << tensor_dim_align << ")" << std::endl
              ;
#endif

    cuda_parallel_launch_local_memory< Multiply ><<< dGrid , dBlock , shmem >>>( *this );
    // Impl::CudaParallelLaunch< Multiply >( *this, dGrid , dBlock , shmem );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT == MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_SHARED_BLOCK_A_X

namespace KokkosArray {
namespace Impl {

/** \brief
 *
 *  finite_elem_row = blockIdx.y ;
 *  stochastic_row  = threadIdx.y + blockDim.y * blockIdx.x ;
 *
 *  Read all of stochastic block coefficients for A and x into shared memory
 *  and perform random accesses into those arrays.
 *
 *  Read all of tensor product for the stochastic row into shared memory.
 */

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensorLegendre< MatrixScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type >  matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >                vector_type ;

private:

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;

public:

  __device__
  void operator()(void) const
  {
    const size_type row_count            = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dim           = m_A.block.dimension();
    const size_type tensor_max_row_width = m_A.block.max_row_width();

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_x[ tensor_dim ]
    //   sh_A[ tensor_dim ]
    //   sh_tvalue[ maximum_tensor_row_width * blockDim.y ]
    //   sh_tcoord[ maximum_tensor_row_width * blockDim.y ]

    volatile VectorScalar * const sh_work = kokkos_impl_cuda_shared_memory<VectorScalar>();
    VectorScalar * const sh_x      = (VectorScalar *)( sh_work + blockDim.x * blockDim.y );
    MatrixScalar * const sh_A      = (MatrixScalar *)( sh_x    + tensor_dim );
    MatrixScalar * const sh_tvalue = (MatrixScalar *)( sh_A    + tensor_dim );
    unsigned     * const sh_tcoord = (unsigned *)( sh_tvalue + tensor_max_row_width * blockDim.y );

    MatrixScalar * const sh_warp_tvalue = sh_tvalue + tensor_max_row_width * threadIdx.y ;
    unsigned     * const sh_warp_tcoord = sh_tcoord + tensor_max_row_width * threadIdx.y ;

    // Block size and thread id within the entire block:

    const size_type nid = CudaTraits::WarpSize * blockDim.y ;
    const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

    // Tensor row for this warp:

    // const size_type iyInner = threadIdx.y + blockDim.y * blockIdx.x ;

    // Each warp's tensor row offsets:

    const unsigned tensor_row_beg   = blockDim.y * blockIdx.x ;
    const unsigned tensor_row_end   = tensor_row_beg + blockDim.y < tensor_dim
                                    ? tensor_row_beg + blockDim.y : tensor_dim ;
    const unsigned tensor_row_count = tensor_row_end - tensor_row_beg ;

    if ( tid < 2 * tensor_row_count + 1 ) {
      sh_tcoord[tid] = m_A.block.m_entry_offset( 2 * tensor_row_beg + tid );
    }
    else if ( tid < 2 * blockDim.y + 1 ) {
      sh_tcoord[tid] = 0 ;
    }

    __syncthreads(); // Wait for read of tensor column offset data

    // Tensor columns for this warp:

    const size_type iBeg       = sh_tcoord[ 2 * threadIdx.y ];
    const size_type iCountDiag = sh_tcoord[ 2 * threadIdx.y + 1 ] > iBeg 
                               ? sh_tcoord[ 2 * threadIdx.y + 1 ] - iBeg : 0 ;
    const size_type iCountAll  = sh_tcoord[ 2 * threadIdx.y + 2 ] > iBeg 
                               ? sh_tcoord[ 2 * threadIdx.y + 2 ] - iBeg : 0 ;

    __syncthreads(); // Wait for query of tensor column offset data

    // Each warp loads its tensor row into shared memory:

    if ( iCountAll ) {
      const MatrixScalar * const a = & m_A.block.m_value(iBeg);
      const unsigned     * const c = & m_A.block.m_coordinate(iBeg);

      for ( size_type i = threadIdx.x ; i < iCountAll ; i += blockDim.x ) {
        sh_warp_tcoord[i] = c[i];
        sh_warp_tvalue[i] = a[i];
      }
    }

    // Each warp computes a row of the stochastic block 
    // for coalesced reads and no need for explicit synchronization.

    for ( size_type iyOuter = blockIdx.y ; iyOuter < row_count ; iyOuter += gridDim.y ) {

      // Loop over columns in the discrete (finite element) system.

      const size_type iBlockEntryEnd = m_A.graph.row_map[ iyOuter + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ iyOuter ];

      VectorScalar y = 0 ;

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        { // Read stochastic coefficients from the vector and matrix into shared memory.
          const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );

          const VectorScalar * const x = & m_x(        0 , iBlockColumn );
          const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

          // Wait for X and A to be used in the previous iteration before reading new values.
          __syncthreads();

          // Coalesced read by the whole block from global memory:
          for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
            sh_x[i] = x[i] ; // m_x(        i , iBlockColumn );
            sh_A[i] = A[i] ; // m_A.values( i , iBlockEntry );
          }

          __syncthreads();
          // Wait for X and A to be read before using these values in the next iteration.
        }

        // Loop through sparse tensor diagonal contributions:

        for ( size_type i = threadIdx.x ; i < iCountDiag ; i += blockDim.x ) {
          const unsigned j = sh_warp_tcoord[i] ;
          y += sh_warp_tvalue[i] * sh_A[j] * sh_x[j] ;
        }

        // Loop through sparse tensor off-diagonal contributions:

        for ( size_type i = threadIdx.x + iCountDiag ; i < iCountAll ; i += blockDim.x ) {
          const unsigned kj = sh_warp_tcoord[i];
          const unsigned j  = kj & 0x0ffff ;
          const unsigned k  = kj >> 16 ;
          y += sh_warp_tvalue[i] * ( sh_A[j] * sh_x[k] + sh_A[k] * sh_x[j] );
        }
      }

      // Reduction of 'y' within the warp

      sh_work[ tid ] = y ;

      if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+16];
      if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 8];
      if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 4];
      if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 2];
      if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 1];

      __syncthreads(); // Wait for warps to complete

      if ( tid < tensor_row_count ) {
        // Coalesced write:
        m_y( tensor_row_beg + tid , iyOuter ) = sh_work[ tid * blockDim.x ];
      }
    }
  }

  //------------------------------------

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
  : m_A( A ), m_x( x ), m_y( y ) {}

  void run() const
  {
    const size_type row_count            = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dim           = m_A.block.dimension();
    const size_type tensor_max_row_width = m_A.block.max_row_width();

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_x[ tensor_dim ]
    //   sh_A[ tensor_dim ]
    //   sh_tvalue[ maximum_tensor_row_width * blockDim.y ]
    //   sh_tcoord[ maximum_tensor_row_width * blockDim.y ]

    const int nsh_warp = sizeof(VectorScalar) * CudaTraits::WarpSize // sh_work
                       + sizeof(MatrixScalar) * tensor_max_row_width // sh_tvalue
                       + sizeof(unsigned)     * tensor_max_row_width // sh_tcoord ;
                       ;

    const int nsh_base = sizeof(VectorScalar) * tensor_dim  // sh_x
                       + sizeof(MatrixScalar) * tensor_dim  // sh_A
                       ;

    enum { ShCapacity = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity };

    int nWarp = ( ShCapacity - nsh_base ) / nsh_warp ;

    if ( nWarp < 1 ) {
      std::ostringstream msg ;
      msg << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply "
          << " FAILED to meet shared memory need: "
          << nsh_base << " + " << nsh_warp
          << " * nWarp <= " << ShCapacity ;
      throw std::runtime_error( msg.str() );
    }

    nWarp = std::min( nWarp , (int) cuda_internal_maximum_warp_count() );
    nWarp = std::min( nWarp , (int) tensor_dim );

    const size_type shmem = nsh_base + nsh_warp * nWarp ;

    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );

    const dim3 dGrid( ( tensor_dim + nWarp - 1 ) / nWarp , row_count , 1 );

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensorLegendre ... > >"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem << ")" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dim(" << tensor_dim << ")" << std::endl
              << "  tensor_max_row_width(" << tensor_max_row_width << ")" << std::endl
              ;
#endif

    cuda_parallel_launch_local_memory< Multiply ><<< dGrid , dBlock , shmem >>>( *this );
    // Impl::CudaParallelLaunch< Multiply >( *this, dGrid , dBlock , shmem );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT == MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_GLOBAL_BLOCK_A_X

namespace KokkosArray {
namespace Impl {

/** \brief
 *
 *  finite_elem_row = blockIdx.y ;
 *  stochastic_row  = threadIdx.y + blockDim.y * blockIdx.x ;
 *
 *  Read all of tensor product for the stochastic row into shared memory.
 *
 *  Read A and x from global memory.
 */

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensorLegendre< MatrixScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type >  matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >                vector_type ;

private:

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;

public:

  __device__
  void operator()(void) const
  {
    const size_type row_count            = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dim           = m_A.block.dimension();
    const size_type tensor_max_row_width = m_A.block.max_row_width();

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_tvalue[ maximum_tensor_row_width * blockDim.y ]
    //   sh_tcoord[ maximum_tensor_row_width * blockDim.y ]

    volatile VectorScalar * const sh_work = kokkos_impl_cuda_shared_memory<VectorScalar>();

    MatrixScalar * const sh_tvalue = (MatrixScalar *)( sh_work + blockDim.x * blockDim.y );
    unsigned     * const sh_tcoord = (unsigned *)( sh_tvalue + tensor_max_row_width * blockDim.y );

    MatrixScalar * const sh_warp_tvalue = sh_tvalue + tensor_max_row_width * threadIdx.y ;
    unsigned     * const sh_warp_tcoord = sh_tcoord + tensor_max_row_width * threadIdx.y ;

    // Block size and thread id within the entire block:

    // const size_type nid = CudaTraits::WarpSize * blockDim.y ;
    const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

    // Tensor row for this warp:

    // const size_type iyInner = threadIdx.y + blockDim.y * blockIdx.x ;

    // Each warp's tensor row offsets:

    const unsigned tensor_row_beg   = blockDim.y * blockIdx.x ;
    const unsigned tensor_row_end   = tensor_row_beg + blockDim.y < tensor_dim
                                    ? tensor_row_beg + blockDim.y : tensor_dim ;
    const unsigned tensor_row_count = tensor_row_end - tensor_row_beg ;

    if ( tid < 2 * tensor_row_count + 1 ) {
      sh_tcoord[tid] = m_A.block.m_entry_offset( 2 * tensor_row_beg + tid );
    }
    else if ( tid < 2 * blockDim.y + 1 ) {
      sh_tcoord[tid] = 0 ;
    }

    __syncthreads(); // Wait for read of tensor column offset data

    // Tensor columns for this warp:

    const size_type iBeg       = sh_tcoord[ 2 * threadIdx.y ];
    const size_type iCountDiag = sh_tcoord[ 2 * threadIdx.y + 1 ] > iBeg 
                               ? sh_tcoord[ 2 * threadIdx.y + 1 ] - iBeg : 0 ;
    const size_type iCountAll  = sh_tcoord[ 2 * threadIdx.y + 2 ] > iBeg 
                               ? sh_tcoord[ 2 * threadIdx.y + 2 ] - iBeg : 0 ;

    __syncthreads(); // Wait for query of tensor column offset data

    // Each warp loads its tensor row into shared memory:

    if ( iCountAll ) {
      const MatrixScalar * const a = & m_A.block.m_value(iBeg);
      const unsigned     * const c = & m_A.block.m_coordinate(iBeg);

      for ( size_type i = threadIdx.x ; i < iCountAll ; i += blockDim.x ) {
        sh_warp_tcoord[i] = c[i];
        sh_warp_tvalue[i] = a[i];
      }
    }

    // Each warp computes a row of the stochastic block 
    // for coalesced reads and no need for explicit synchronization.

    for ( size_type iyOuter = blockIdx.y ; iyOuter < row_count ; iyOuter += gridDim.y ) {

      // Loop over columns in the discrete (finite element) system.

      const size_type iBlockEntryEnd = m_A.graph.row_map[ iyOuter + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ iyOuter ];

      VectorScalar y = 0 ;

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        // Read stochastic coefficients from the vector and matrix into shared memory.
        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );

        const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

        // Loop through sparse tensor diagonal contributions:

        for ( size_type i = threadIdx.x ; i < iCountDiag ; i += blockDim.x ) {
          const unsigned j = sh_warp_tcoord[i] ;
          y += sh_warp_tvalue[i] * A[j] * x[j] ;
        }

        // Loop through sparse tensor off-diagonal contributions:

        for ( size_type i = threadIdx.x + iCountDiag ; i < iCountAll ; i += blockDim.x ) {
          const unsigned kj = sh_warp_tcoord[i];
          const unsigned j  = kj & 0x0ffff ;
          const unsigned k  = kj >> 16 ;
          y += sh_warp_tvalue[i] * ( A[j] * x[k] + A[k] * x[j] );
        }
      }

      // Reduction of 'y' within the warp

      sh_work[ tid ] = y ;

      if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+16];
      if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 8];
      if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 4];
      if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 2];
      if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_work[tid] += sh_work[tid+ 1];

      __syncthreads(); // Wait for warps to complete

      if ( tid < tensor_row_count ) {
        // Coalesced write:
        m_y( tensor_row_beg + tid , iyOuter ) = sh_work[ tid * blockDim.x ];
      }
    }
  }

  //------------------------------------

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
  : m_A( A ), m_x( x ), m_y( y ) {}

  void run() const
  {
    const size_type row_count            = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dim           = m_A.block.dimension();
    const size_type tensor_max_row_width = m_A.block.max_row_width();

    // Shared memory:
    //   sh_work[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_tvalue[ maximum_tensor_row_width * blockDim.y ]
    //   sh_tcoord[ maximum_tensor_row_width * blockDim.y ]

    const int nsh_warp = sizeof(VectorScalar) * CudaTraits::WarpSize // sh_work
                       + sizeof(MatrixScalar) * tensor_max_row_width // sh_tvalue
                       + sizeof(unsigned)     * tensor_max_row_width // sh_tcoord ;
                       ;

    const int nsh_base = 0 ;

    enum { ShCapacity = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity };

    int nWarp = ( ShCapacity - nsh_base ) / nsh_warp ;

    if ( nWarp < 1 ) {
      std::ostringstream msg ;
      msg << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply "
          << " FAILED to meet shared memory need: "
          << nsh_base << " + " << nsh_warp
          << " * nWarp <= " << ShCapacity ;
      throw std::runtime_error( msg.str() );
    }

    nWarp = std::min( nWarp , (int) cuda_internal_maximum_warp_count() );
    nWarp = std::min( nWarp , (int) tensor_dim );

    const size_type shmem = nsh_base + nsh_warp * nWarp ;

    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );

    const dim3 dGrid( ( tensor_dim + nWarp - 1 ) / nWarp , row_count , 1 );

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensorLegendre ... > >"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem << ")" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dim(" << tensor_dim << ")" << std::endl
              << "  tensor_max_row_width(" << tensor_max_row_width << ")" << std::endl
              ;
#endif

    cuda_parallel_launch_local_memory< Multiply ><<< dGrid , dBlock , shmem >>>( *this );
    // Impl::CudaParallelLaunch< Multiply >( *this, dGrid , dBlock , shmem );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT == MULTIPLY_CRS_PRODUCT_TENSOR_LEGENDRE_VARIANT_VARIABLE_BLOCK_A_X

namespace KokkosArray {
namespace Impl {

/** \brief
 *
 *  finite_elem_row = blockIdx.y ;
 *  stochastic_row  =  + gang_size[ blockIdx.x ]
 *
 *  Read all of tensor product for the stochastic row into shared memory.
 *
 *  Read A and x from global memory.
 */

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensorLegendre< MatrixScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type >  matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >                vector_type ;

private:

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;

  View< unsigned * , Cuda > m_tensor_row_offset ;
  size_type m_shmem ;

  enum { WarpCount   = 8 };
  enum { ThreadCount = WarpCount * CudaTraits::WarpSize };

public:

  __device__
  void operator()(void) const
  {
    // Finite element and stochastic dimensions:
    const size_type fem_row_count  = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_dim     = m_A.block.dimension();

    // This block's span of tensor rows:
    const size_type tensor_row_begin = m_tensor_row_offset( blockIdx.x );
    const size_type tensor_row_next  = m_tensor_row_offset( blockIdx.x + 1 );
    const size_type tensor_row_end   = tensor_row_next < tensor_dim
                                     ? tensor_row_next : tensor_dim ;

    // All threads within this block have the same gang_size
    // which is a multiple of 2 and gang_size <= WarpSize.

    const size_type gang_count = tensor_row_next - tensor_row_begin ;
    const size_type gang_size  = blockDim.x  / gang_count ;
    const size_type gang_rank  = threadIdx.x / gang_size ;
    const size_type gang_tid   = threadIdx.x % gang_size ;

    // Shared memory:
    //   sh_work[ blockDim.x ]
    //   sh_tvalue[ max_entry_count ]
    //   sh_tcoord[ max_entry_count ]

    volatile VectorScalar * const sh_work = kokkos_impl_cuda_shared_memory<VectorScalar>();
    unsigned * const sh_toffset = (unsigned *)( sh_work );

    // Read the offset array from the crs tensor:

    if ( threadIdx.x <= 2 * gang_count ) {
      const size_type jEnd = 2 * tensor_row_end ;
      const size_type j    = 2 * tensor_row_begin + threadIdx.x ;
      sh_toffset[ threadIdx.x ] = m_A.block.m_entry_offset( j < jEnd ? j : jEnd );
    }

    __syncthreads();

    const size_type tensor_entry_begin = sh_toffset[0] ;
    const size_type tensor_entry_count = sh_toffset[ 2 * gang_count ] - tensor_entry_begin ;
    const size_type tensor_entry_gang_begin = sh_toffset[ 2 * gang_rank ] ;
    const size_type tensor_entry_gang_count = sh_toffset[ 2 * ( gang_rank + 1 ) ] - tensor_entry_gang_begin ;

    const size_type gang_limit = gang_size < tensor_entry_gang_count 
                               ? gang_size : tensor_entry_gang_count ;

    MatrixScalar * const sh_tvalue  = (MatrixScalar *)( sh_work + blockDim.x );
    unsigned     * const sh_tcoord  = (unsigned *)( sh_tvalue + tensor_entry_count );

    MatrixScalar * const sh_gang_tvalue = sh_tvalue + tensor_entry_gang_begin - tensor_entry_begin ;
    unsigned     * const sh_gang_tcoord = sh_tcoord + tensor_entry_gang_begin - tensor_entry_begin ;

    // Read coordinates and values from crs tensor:
    {
      const MatrixScalar * const a = & m_A.block.m_value( tensor_entry_begin );
      const unsigned     * const c = & m_A.block.m_coordinate( tensor_entry_begin );

      for ( size_type i = threadIdx.x ; i < tensor_entry_count ; i += blockDim.x ) {
        sh_tvalue[i] = a[i] ;
        sh_tcoord[i] = c[i] ;
      }
    }

    __syncthreads(); // Wait for read of tensor data

    // Each gang computes a row of the stochastic block 
    // for coalesced reads and no need for explicit synchronization.

    for ( size_type iyOuter = blockIdx.y ; iyOuter < fem_row_count ; iyOuter += gridDim.y ) {

      // Loop over columns in the discrete (finite element) system.

      const size_type iBlockEntryEnd = m_A.graph.row_map[ iyOuter + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ iyOuter ];

      VectorScalar y = 0 ;

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        // Read stochastic coefficients from the vector and matrix into shared memory.
        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );

        const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

        // Loop through sparse tensor off-diagonal contributions:

        for ( size_type i = gang_tid ; i < tensor_entry_gang_count ; i += gang_size ) {
          const unsigned kj = sh_gang_tcoord[i];
          const unsigned j  = kj & 0x0ffff ;
          const unsigned k  = kj >> 16 ;
          y += sh_gang_tvalue[i] * ( A[j] * x[k] + A[k] * x[j] );
        }
      }

      // Reduction of 'y' within the gang (power of two and <= WarpSize)

      sh_work[ threadIdx.x ] = y ;

      if ( gang_tid +  1 < gang_limit ) {
      if ( gang_tid +  2 < gang_limit ) {
      if ( gang_tid +  4 < gang_limit ) {
      if ( gang_tid +  8 < gang_limit ) {
      if ( gang_tid + 16 < gang_limit ) {
        sh_work[ threadIdx.x ] += sh_work[ threadIdx.x + 16 ]; }
        sh_work[ threadIdx.x ] += sh_work[ threadIdx.x +  8 ]; }
        sh_work[ threadIdx.x ] += sh_work[ threadIdx.x +  4 ]; }
        sh_work[ threadIdx.x ] += sh_work[ threadIdx.x +  2 ]; }
        sh_work[ threadIdx.x ] += sh_work[ threadIdx.x +  1 ]; }

      __syncthreads(); // Wait for warps to complete

      if ( tensor_row_begin + threadIdx.x < tensor_row_end ) {
        const size_type b = threadIdx.x * gang_size ;

        // If the gang size exceeds the warp size must reduce across warps.
        for ( size_type n = 32 ; n < gang_limit ; n += 32 ) {
          sh_work[ b ] += sh_work[ b + n ];
        }

        // Coalesced write:
        m_y( tensor_row_begin + threadIdx.x , iyOuter ) = sh_work[ b ];
      }
    }
  }

  //------------------------------------

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
  : m_A( A ), m_x( x ), m_y( y )
  , m_tensor_row_offset()
  , m_shmem(0)
  {
    enum { ShCapacity = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity };

    const size_type target_shmem   = ShCapacity / 4 ;
    const size_type tensor_dim     = m_A.block.dimension();
    const size_type max_gang_count = ThreadCount / 2 ;

    // Analyze the tensor to partition rows into blocks of roughly equal
    // nonzero count.
    // Limit the size of a gang to a warp for fast reductions.

    typename tensor_type::array_unsigned_host_type
      tensor_entry_offset = create_mirror( m_A.block.m_entry_offset );

    deep_copy( tensor_entry_offset , m_A.block.m_entry_offset );

    std::vector<unsigned> tensor_row_offset(1,0u);

    for ( size_type tensor_row = 0 ; tensor_row < tensor_dim ; ) {
      // Power of two gang_size thus a power of two gang_count
      // so that gang_size * gang_count == ThreadCount

      const size_type row_entry_begin = tensor_entry_offset(2*tensor_row);

      size_type gang_count = 1 ;

      for ( ; gang_count < max_gang_count &&
              tensor_row + gang_count < tensor_dim ; gang_count <<= 1 ) {

        const size_type tensor_row_end = std::min( tensor_dim , tensor_row + ( gang_count << 1 ) );

        const size_type row_entry_end   = tensor_entry_offset(2*tensor_row_end);
        const size_type row_entry_count = row_entry_end - row_entry_begin ;
        const size_type block_shmem
          = sizeof(VectorScalar) * ThreadCount    // sh_work[ threads ]
          + sizeof(MatrixScalar) * row_entry_count // sh_tvalue[ entries ]
          + sizeof(unsigned)     * row_entry_count // sh_tcoord[ entries ]
          ;

        if ( target_shmem < block_shmem ) break ;
      }

      const size_type tensor_row_next = tensor_row + gang_count ;

      {
        const size_type tensor_row_end  = std::min( tensor_dim , tensor_row_next );
        const size_type row_entry_end   = tensor_entry_offset(2*tensor_row_end);
        const size_type row_entry_count = row_entry_end - row_entry_begin ;

        const size_type block_shmem
          = sizeof(VectorScalar) * ThreadCount    // sh_work[ threads ]
          + sizeof(MatrixScalar) * row_entry_count // sh_tvalue[ entries ]
          + sizeof(unsigned)     * row_entry_count // sh_tcoord[ entries ]
          ;

        m_shmem = std::max( m_shmem , block_shmem );
      }

      tensor_row = tensor_row_next ;

      tensor_row_offset.push_back( tensor_row );
    }

    m_tensor_row_offset = View< unsigned * , Cuda >( "tensor_row_work_offset" , tensor_row_offset.size() );

    View< unsigned * , Cuda >::HostMirror h_tensor_row_offset = create_mirror( m_tensor_row_offset );

    for ( size_type i = 0 ; i < h_tensor_row_offset.dimension_0() ; ++i ) {
      h_tensor_row_offset(i) = tensor_row_offset[i] ;
    }

    deep_copy( m_tensor_row_offset , h_tensor_row_offset );

#if 0

    const size_type fem_row_count  = m_A.graph.row_map.dimension_0() - 1 ;

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensorLegendre ... > >"
              << std::endl 
              << "  grid(" << m_tensor_row_offset.dimension_0() - 1
              << "," << fem_row_count << ")" << std::endl
              << "  block(" << ThreadCount << ")" << std::endl
              << "  shmem(" << m_shmem << ")" << std::endl
              << "  fem_row_count(" << fem_row_count << ")" << std::endl
              << "  tensor_dim(" << tensor_dim << ")" << std::endl
              << "  tensor_blocking("
              ;
    for ( unsigned i = 0 ; i < tensor_row_offset.size() ; ++i ) {
      std::cout << " ( " << tensor_row_offset[i]
                << " : " << tensor_entry_offset[i*2] 
                << " )" ;
    }
    std::cout << " )" << std::endl ;
#endif

  }

  void run() const
  {
    const size_type fem_row_count  = m_A.graph.row_map.dimension_0() - 1 ;
    const size_type tensor_subsets = m_tensor_row_offset.dimension_0() - 1 ;

    const dim3 dBlock( ThreadCount , 1 , 1 );
    const dim3 dGrid( tensor_subsets , fem_row_count , 1 );

    cuda_parallel_launch_local_memory< Multiply ><<< dGrid , dBlock , m_shmem >>>( *this );
    // Impl::CudaParallelLaunch< Multiply >( *this, dGrid , dBlock , shmem );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_CUDA_CRSPRODUCTTENSORLEGENDRE_HPP */

