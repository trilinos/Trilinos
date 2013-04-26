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

#ifndef KOKKOSARRAY_CUDA_SPARSEPRODUCTTENSORLEGENDRE_HPP
#define KOKKOSARRAY_CUDA_SPARSEPRODUCTTENSORLEGENDRE_HPP

#include <utility>
#include <sstream>
#include <stdexcept>
#include <Kokkos_Atomic.hpp>
#include <KokkosArray_SparseProductTensorLegendre.hpp>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>


namespace KokkosArray {
namespace Impl {

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix<
    SparseProductTensorLegendre< TensorScalar , Cuda ,
                                 SparseProductTensorLegendreVariant_Default > ,
    MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda device_type ;
  typedef typename device_type::size_type  size_type ;

private:


  typedef SparseProductTensorLegendre< MatrixScalar , device_type ,
                                       SparseProductTensorLegendreVariant_Default > 
    tensor_type ;

  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type >  matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , device_type >         vector_type ;

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;
  int   m_column_count ;
  int   m_shared_size ;

public:

  Multiply( const matrix_type & arg_A ,
            const vector_type & arg_x ,
            const vector_type & arg_y )
  : m_A( arg_A )
  , m_x( arg_x )
  , m_y( arg_y )
  , m_column_count(0)
  , m_shared_size(0)
  {
    // Use at most a quarter of shared memory to allow 4 blocks per SMP
    enum { SharedMemoryLimit = Impl::CudaTraits::SharedMemoryCapacity / 4 };

    const size_type shared_size_per_column =
      m_A.block.dimension() * ( sizeof(VectorScalar) + sizeof(MatrixScalar) );

    m_column_count = SharedMemoryLimit / shared_size_per_column ;
    m_shared_size  = CudaTraits::warp_align( m_column_count * m_A.block.dimension() ) *
                     ( sizeof(VectorScalar) + sizeof(MatrixScalar) );

    while ( SharedMemoryLimit < m_shared_size ) {
      m_column_count-- ;
      m_shared_size  = CudaTraits::warp_align( m_column_count * m_A.block.dimension() ) *
                       ( sizeof(VectorScalar) + sizeof(MatrixScalar) );
    }

    if ( m_A.block.dimension() != m_x.dimension_0() ||
         m_A.block.dimension() != m_y.dimension_0() ||
         m_A.block.dimension() != m_A.values.dimension_0() ||
         m_x.dimension_1() != m_y.dimension_1() ||
         0 == m_column_count ) {
      throw std::runtime_error("KokkosArray::Multiply< BlockCrsMatrix< SparseProductTensor<> ... > > dimension error");
    }
  }

  //  Maximum number of tensor entries that are processed by a thread.
  //  The number of registers consumed by the thread is proportional.
  enum { ThreadMaxTensorEntry = 4 };

  void run() const
  {
    // Block and warp residency is impacted by the maximum threads per block
    enum { MAX_THREAD = 8 * Impl::CudaTraits::WarpSize };
    enum { MAX_TENSOR = Impl::CudaTraits::UpperBoundGridCount * ThreadMaxTensorEntry * MAX_THREAD };

    if ( MAX_TENSOR < m_A.block.entry_count() ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Multiply< BlockCrsMatrix< SparseProductTensor<> ... > > "
          << " Tensor too large for Cuda implementation" ;
      msg << " : tensor-entry-count = " << m_A.block.entry_count() << " > " << MAX_TENSOR ;
      throw std::runtime_error( msg.str() );
    }

    // Maximum number of tensor entries that can be managed by a thread block.
    //   Need enough blocks in the 'grid.x' dimension to hold all tensor entries.
    //   There will be 'grid.x' redundant reads of A and X values.
    const size_type maxTensorBlock = std::min( unsigned( ThreadMaxTensorEntry * MAX_THREAD ) ,
                                               unsigned( m_A.block.entry_count() ) );

    // There are 'grid.y' redundant reads of the tensor.
    // However, this is the primary source of block-level parallelism.

    const dim3 grid( ( ( m_A.block.entry_count() + maxTensorBlock - 1 ) / maxTensorBlock ) ,
                     std::min( unsigned(m_y.dimension_1()) , unsigned(Impl::CudaTraits::UpperBoundGridCount) ) , 1 );

    const dim3 block( MAX_THREAD / 2 , 2 , 1 );

    //------------------------------------------------------------------------

#if 0

    std::cout << "KokkosArray::Multiply< BlockCrsMatrix< SparseProductTensorLegendre< ... > > >::run() "
              << std::endl
              << "  tensor_dimen = " << m_A.block.dimension()   << std::endl
              << "  tensor_entry = " << m_A.block.entry_count() << std::endl
              << "  column_count = " << m_column_count << std::endl
              << "  grid  = " << grid.x  << " , " << grid.y << std::endl
              << "  block = " << block.x << " , " << block.y << std::endl
              << "  shmem = " << m_shared_size << std::endl ;

#endif

    //------------------------------------------------------------------------

    Impl::cuda_parallel_launch_local_memory<<< grid , block , m_shared_size >>>( *this );
  }

  // This thread is reponsible for up to ThreadMaxTensorEntry tensor entries.
  //
  // Threads apply multiple tensor entries to reduce the amount of redundant
  // reads of the matrix and vector values.
  //
  // The amount of redundancy is gridDim.x where gridDim.x is large enough
  // to satisfy:
  //   m_A.block.entry_count() <= ( ThreadMaxTensorEntry * blockDim.x * blockDim.y * gridDim.x )
  //
  //
  // The tensor is redundantly read by thread blocks in the gridDim.y axis.
  // These blocks are associated with the outer dimension of the multivector.
  // Both parallelism and redundant tensor reads increase with gridDim.y;
  // thus gridDim.y should be chosen to maximize utilization of the hardware
  // but no larger than needed to meet this condition.

  __device__
  void operator()(void) const
  {
    const int tensor_dim = m_A.block.dimension();
    const int sTensorEntry = blockDim.x * blockDim.y * gridDim.x ; /* Stride for tensor entries */
    const int iTensorEntry = threadIdx.x + blockDim.x * ( threadIdx.y + blockDim.y * blockIdx.x );

    VectorScalar * const shX = kokkos_impl_cuda_shared_memory<VectorScalar>();
    MatrixScalar * const shA = (MatrixScalar *) ( shX + CudaTraits::warp_align( m_column_count * tensor_dim ) );

    VectorScalar yi0[ ThreadMaxTensorEntry ] , yi1[ ThreadMaxTensorEntry ] ;
    VectorScalar yj0[ ThreadMaxTensorEntry ] , yj1[ ThreadMaxTensorEntry ] ;
    VectorScalar yk0[ ThreadMaxTensorEntry ] , yk1[ ThreadMaxTensorEntry ] ;

    int i[ ThreadMaxTensorEntry ] ;
    int j[ ThreadMaxTensorEntry ] ;
    int k[ ThreadMaxTensorEntry ] ;
    int n_active = 0 ;

    //------------------------------------
    // Read tensor coordinates which this thread will process.

    for ( int iTensor = iTensorEntry ;
              iTensor < m_A.block.entry_count() && n_active < ThreadMaxTensorEntry ;
              iTensor += sTensorEntry , ++n_active ) {

      const unsigned long ic = m_A.block.m_coordinate( iTensor );

      i[n_active] = ( ic >> 32 ) & 0x0ffff ;
      j[n_active] = ( ic >> 16 ) & 0x0ffff ;
      k[n_active] = ( ic       ) & 0x0ffff ;
    }

    //------------------------------------
    // Iterate through outer rows of the linear system.

    for ( int iBlockRow = blockIdx.y ;
              iBlockRow < m_y.dimension_1() ;
              iBlockRow += gridDim.y ) {

      for ( int n = 0 ; n < n_active ; ++n ) {
        yi0[n] = 0 ; yi1[n] = 0 ;
        yj0[n] = 0 ; yj1[n] = 0 ;
        yk0[n] = 0 ; yk1[n] = 0 ;
      }

      // Iterate through outer columns of the linear system:

      const size_type iBlockEntryEnd = m_A.graph.row_map[iBlockRow+1];
            size_type iBlockEntry    = m_A.graph.row_map[iBlockRow];

      for ( ; iBlockEntry < iBlockEntryEnd ; iBlockEntry += m_column_count ) {

        // Load a span of inner stochastic coefficients into shared memory:

        const size_type total_coeff = tensor_dim *
           ( iBlockEntry + m_column_count < iBlockEntryEnd ? m_column_count : iBlockEntryEnd - iBlockEntry );

        __syncthreads(); // Wait for previous load to be used.

        // Half of the threads read A, other half read X
        // Maximize thread utilization for the global memory read.
        if ( 0 == threadIdx.y ) {
          for ( int n = threadIdx.x ; n < total_coeff ; n += blockDim.x ) {
            shA[n] = m_A.values( n % tensor_dim , iBlockEntry + n / tensor_dim );
          }
        }
        else {
          for ( int n = threadIdx.x ; n < total_coeff ; n += blockDim.x ) {
            const size_type iBlockColumn = m_A.graph.entries( iBlockEntry + n / tensor_dim );
            shX[n] = m_x( n % tensor_dim , iBlockColumn );
          }
        }

        __syncthreads(); // Wait for current load to complete.

        // Local stochastic block multiplication:
        // Faster to loop over active tensor entries and then outer columns.

        for ( int n = 0 ; n < n_active ; ++n ) {
          for ( int in = i[n] , jn = j[n] , kn = k[n] ; in < total_coeff ;
                in += tensor_dim ,
                jn += tensor_dim ,
                kn += tensor_dim ) {

            yi0[n] += shA[ jn ] * shX[ kn ] ; yi1[n] += shX[ jn ] * shA[ kn ] ;
            yj0[n] += shA[ kn ] * shX[ in ] ; yj1[n] += shX[ kn ] * shA[ in ] ;
            yk0[n] += shA[ in ] * shX[ jn ] ; yk1[n] += shX[ in ] * shA[ jn ] ;
          }
        }
      }

      // Add this thread's contribution to the multiplication:
      // Guaranteed by construction that: i >= j >= k
      // When: i != j  then  i >  j >= k  and therefore  i != k
      // When: j != k  then  i >= j >  k  and therefore  i != k

      // Atomic count = nTensorEntry * 3 , and very imbalanced...

      for ( int n = 0 , iTensor = iTensorEntry ; n < n_active ; ++n , iTensor += sTensorEntry ) {

        const TensorScalar v = m_A.block.m_value( iTensor );
        const int neq_ij = i[n] != j[n] ;
        const int neq_jk = j[n] != k[n] ;

                      Kokkos::atomic_fetch_add( & m_y( i[n] , iBlockRow ) , v * ( yi0[n] + yi1[n] * neq_jk ) );
        if ( neq_ij ) Kokkos::atomic_fetch_add( & m_y( j[n] , iBlockRow ) , v * ( yj0[n] + yj1[n] ) );
        if ( neq_jk ) Kokkos::atomic_fetch_add( & m_y( k[n] , iBlockRow ) , v * ( yk0[n] + yk1[n] * neq_ij ) );
      }
    }
  }
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_CUDA_SPARSEPRODUCTTENSORLEGENDRE_HPP */

