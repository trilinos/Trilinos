/*
//@HEADER
// ************************************************************************
// 
//                         Kokkos Array
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

#ifndef KOKKOS_CUDA_PRODUCTTENSOR_HPP
#define KOKKOS_CUDA_PRODUCTTENSOR_HPP

#include <utility>
#include <sstream>
#include <stdexcept>
#include <Cuda/Kokkos_Cuda_Parallel.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< typename TensorScalar >
class Multiply< SparseProductTensor< 3 , TensorScalar , Cuda > , void , void >
{
public:
  typedef Cuda::size_type size_type ;
  typedef SparseProductTensor< 3 , TensorScalar , Cuda > block_type ;

  static size_type matrix_size( const block_type & tensor )
  { return tensor.dimension(); }
  
  static size_type vector_size( const block_type & tensor )
  { return tensor.dimension(); }

  //----------------------------------------

  static dim3 thread_block( const block_type & block )
  {
    const int n = cuda_internal_maximum_warp_count() * CudaTraits::WarpSize ;

    if ( n < block.dimension() ) {
      throw std::runtime_error( std::string("Kokkos::Impl::Multiply< SparseProductTensor<3,Scalar,Cuda> > ERROR: block too large") );
    }

    const dim3 d( block.dimension() ,
                  std::min( n / block.dimension() , block.entry_count() ) , 1 );

    return d ;
  }

  template< typename VectorValue >
  __host__
  static size_type shmem_size( const block_type & block )
  {
    const dim3 d = thread_block( block );

    const size_type size_reduce = d.x * d.y * sizeof(VectorValue);

    const size_type size_work =
      2 * d.x * sizeof(VectorValue) +
      d.x * d.y * sizeof(TensorScalar) +
      d.x * d.y * 3 * sizeof(int);

    return std::max( size_reduce , size_work );
  }

  //----------------------------------------

  template< typename MatrixValue , typename VectorValue >
  __device__
  static VectorValue apply( const block_type & block ,
                            const MatrixValue * const a ,
                            const VectorValue * const x )
  {
    VectorValue  * const shX = kokkos_impl_cuda_shared_memory<VectorValue>();
    VectorValue  * const shA = (VectorValue  *)( shX + blockDim.x );

    VectorValue yVal = 0 ;

    // Load shared x and A vectors:

    if ( threadIdx.y == 0 ) {
      shX[ threadIdx.x ] = x[ threadIdx.x ];
    }

    if ( threadIdx.y == blockDim.y - 1 ) {
      shA[ threadIdx.x ] = a[ threadIdx.x ];
    }

    __syncthreads(); // wait for all data reads

    // Loop through sparse tensor contributions:

#if 0

    const size_type nThread = blockDim.x * blockDim.y ;
    const size_type iThread = threadIdx.x + blockDim.x * threadIdx.y ;

    TensorScalar * const shValue = (TensorScalar *)( shA + blockDim.x );
    int          * const shCoord = (int *)( shValue + nThread );

    size_type nOuterLoop = ( block.entry_count() + nThread - 1 ) / nThread ;

    for ( size_type iOuterLoop = 0 ;
                    iOuterLoop < nOuterLoop ; ++iOuterLoop ) {

      const size_type base = iOuterLoop * nThread ;
      const size_type iThreadBase = iThread + base ;

      shCoord[ iThread ] = -1 ;

      if ( iThreadBase < block.entry_count() ) {
        shCoord[ iThread               ] = block.coord( iThreadBase , 0 );
        shCoord[ iThread + nThread     ] = block.coord( iThreadBase , 1 );
        shCoord[ iThread + nThread * 2 ] = block.coord( iThreadBase , 2 );
        shValue[ iThread               ] = block.value( iThreadBase );
      }

      __syncthreads(); // Wait for load of 'nThread' entries.

      size_type endEntry = base + nThread ;
      if ( block.entry_count() < endEntry ) endEntry = block.entry_count();

      for ( size_type iEntry = threadIdx.y ;
                      iEntry < endEntry ; iEntry += blockDim.y ) {

        int i = shCoord[ iEntry ];
        int j = shCoord[ iEntry + nThread ];
        int k = shCoord[ iEntry + nThread * 2 ];

        if      ( threadIdx.x == i ) { ; }
        else if ( threadIdx.x == j ) { j = i ; }
        else if ( threadIdx.x == k ) { k = i ; }
        else { i = -1 ; }

        if ( i != -1 ) {
          VectorValue tmp = shX[j] * shA[k];
          if ( j != k ) { tmp += shX[k] * shA[j]; }
          yVal += tmp * shValue[ iEntry ];
        }
      }

      __syncthreads(); // Wait for use of 'nThread' entries.

    }

#else

    for ( size_type iEntry = threadIdx.y ;
                    iEntry < block.entry_count() ; iEntry += blockDim.y ) {

      int i = block.coord(iEntry,0);
      int j = block.coord(iEntry,1);
      int k = block.coord(iEntry,2);

      if      ( threadIdx.x == i ) { ; }
      else if ( threadIdx.x == j ) { j = i ; }
      else if ( threadIdx.x == k ) { k = i ; }
      else { i = -1 ; }

      if ( i != -1 ) {
        VectorValue tmp = shX[j] * shA[k];
        if ( j != k ) { tmp += shX[k] * shA[j]; }
        yVal += tmp * block.value(iEntry);
      }
    }

#endif

    if ( 1 < blockDim.y ) {
      __syncthreads(); // wait for all data uses

      if ( 0 < threadIdx.y ) {
        shX[ threadIdx.x + blockDim.x * threadIdx.y ] = yVal ;
      }

      __syncthreads(); // wait for all data loads

      for ( int i = 1 ; i < blockDim.y ; ++i ) {
        yVal += shX[ threadIdx.x + i * blockDim.x ];
      }
    }

    return yVal ;
  }
};

//----------------------------------------------------------------------------

#if 1

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< SparseProductTensor< 3 , TensorScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  MultiVector< VectorScalar , Cuda > ,
  MultiVector< VectorScalar , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef SparseProductTensor< 3 , TensorScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef MultiVector< VectorScalar , device_type>                 vector_type ;

  enum { WARP_SIZE = Impl::CudaTraits::WarpSize };

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;
  const size_type   m_shared_vector_count ;

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
    : m_A( A ), m_x( x ), m_y( y )
    , m_shared_vector_count( CudaTraits::warp_align( A.block.dimension() ) )
    {}

  inline
  __device__
  static void swap( size_type & i , size_type & j )
  { const size_type tmp = i ; i = j ; j = tmp ; }


  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    // Six warps required for the six uses of a product tensor value
    // 'N' warps required for reading a block of values.
    const size_type nWarp =
      std::max( size_type(6) , Impl::CudaTraits::warp_count( A.block.dimension() ) );

    const size_type nGridX =
      Impl::CudaTraits::warp_count( A.block.entry_count() );

    const size_type shared_size =
      std::max( size_type( WARP_SIZE * ( 3 * sizeof(size_type) + sizeof(VectorScalar) ) ),
                size_type( CudaTraits::warp_align( A.block.dimension() ) * 2 * sizeof(VectorScalar) ) );

    //------------------------------------------------------------------------
    // Sizing error conditions

    if ( nWarp  > Impl::cuda_internal_maximum_warp_count() ||
         nGridX > Impl::CudaTraits::UpperBoundGridCount ||
         shared_size >
            Impl::cuda_internal_maximum_shared_words() * sizeof(size_type) ) {

      std::ostringstream msg ;
      msg << "Kokkos::Multiply< BlockCrsMatrix< SparseProductTensor<> ... > > "
          << " Tensor too large for Cuda implementation" ;
      if ( nWarp  > Impl::cuda_internal_maximum_warp_count() ) {
        msg << " : tensor-dimension = " << A.block.dimension()
            << " > " << Impl::cuda_internal_maximum_warp_count() * WARP_SIZE ;
      }
      if ( nGridX > Impl::CudaTraits::UpperBoundGridCount ) {
        msg << " : tensor-entry-count = " << A.block.entry_count()
            << " > " << Impl::CudaTraits::UpperBoundGridCount * WARP_SIZE ;
      }
      if ( shared_size >
             Impl::cuda_internal_maximum_shared_words() * sizeof(size_type) ) {
        msg << " : requested shared_memory " << shared_size << " > " 
            <<  Impl::cuda_internal_maximum_shared_words() * sizeof(size_type);
      }
      throw std::runtime_error( msg.str() );
    }

    //------------------------------------------------------------------------

    const size_type nGridY =
       std::min( A.graph.row_count() ,
                 Impl::CudaTraits::UpperBoundGridCount / nGridX );

    const dim3 block( WARP_SIZE , nWarp , 1 );
    const dim3 grid(  nGridX , nGridY , 1 );

    Impl::cuda_parallel_launch_local_memory<<< grid , block , shared_size >>>
      ( Multiply( A , x , y ) );
  }

  inline
  __device__
  size_type & shared_tensor_coord( size_type C ) const
  {
    return kokkos_impl_cuda_shared_memory<size_type>()
           [ threadIdx.x + C * WARP_SIZE ];
  }

  inline
  __device__
  VectorScalar & shared_tensor_value() const
  {
    enum { offset = (3 * WARP_SIZE * sizeof(size_type)) / sizeof(VectorScalar) };
    return kokkos_impl_cuda_shared_memory<VectorScalar>()
           [ threadIdx.x + offset ];
  }

  inline
  __device__
  VectorScalar & shared_matrix_value( size_type i ) const
  {
    return kokkos_impl_cuda_shared_memory<VectorScalar>()[i];
  }

  inline
  __device__
  VectorScalar & shared_vector_value( size_type j ) const
  {
    return kokkos_impl_cuda_shared_memory<VectorScalar>()
             [m_shared_vector_count+j];
  }

  // Required:
  //   m_A.block.dimension() == m_A.values.dimension(0)
  //   blockDim.x == WARP_SIZE
  //   blockDim.y >= 6
  //   blockDim.x * blockDim.y >= m_A.block.dimension()
  //   blockDim.x * gridDim.x  >= m_A.block.entry_count();
  //
  //   shared_tensor_size =
  //     WARP_SIZE * ( 3 sizeof(size_type) + sizeof(VectorScalar) )
  //
  //   shared_coeff_size =
  //     m_A.block.dimension() * 2 * sizeof(VectorScalar)
  //                        
  //   shmem_size >= max( shared_tensor_size , shared_coeff_size )
  //
  __device__
  void execute_on_device() const
  {
    // Coalesced read of A and x coefficients into shared memory:
    const size_type iCoeff = threadIdx.x + WARP_SIZE * threadIdx.y ;
    const bool read_coeff  = iCoeff < m_A.block.dimension();

    //------------------------------------------------------------
    // This thread is responsible for one sparse tensor entry

    bool have_tensor = false ;

    {
      // This thread's tensor entry to read into shared memory:
      const size_type iTensor = threadIdx.x + WARP_SIZE * blockIdx.x ;

      if ( ( have_tensor = iTensor < m_A.block.entry_count() ) ) {
        if ( threadIdx.y < 3 ) {
          shared_tensor_coord( threadIdx.y ) =
            m_A.block.coord( iTensor , threadIdx.y );
        }
        else if ( threadIdx.y == 3 ) {
          shared_tensor_value() = m_A.block.value( iTensor );
        }
      }
    }

    __syncthreads();

    const VectorScalar tval = shared_tensor_value();

    // Product tensor symmetry gives six potential uses of
    // the tensor value.  Have six warps, each warp is responsible
    // for applying a non-redundant use of a tensor value.
    //
    //   y(kY,row) += tensor(iA,jX,kY) * A(iA,row,col) * x(jX,col)
    //
    size_type iA = shared_tensor_coord(0);
    size_type jX = shared_tensor_coord(1);
    size_type kY = shared_tensor_coord(2);

    if ( have_tensor ) {
      switch( threadIdx.y ) {
      case 0 :
        // y(k,row) += tensor(i,j,k) * ( A(i,row,col) * x(j,col) )
        have_tensor = true ;
        break ;
      case 1 :
        // y(k,row) += tensor(i,j,k) * ( A(j,row,col) * x(i,col) )
        if ( ( have_tensor = iA != jX ) ) { swap( jX , iA ); }
        break ;
      case 2 :
        // y(j,row) += tensor(i,j,k) * ( A(i,row,col) * x(k,col) )
        if ( ( have_tensor = jX != kY ) ) { swap( kY , jX ); }
        break ;
      case 3 :
        // y(j,row) += tensor(i,j,k) * ( A(k,row,col) * x(i,col) )
        if (( have_tensor = jX != kY && kY != iA )) {
          swap(kY,jX); swap(jX,iA);
        }
        break ;
      case 4 :
        // y(i,row) += tensor(i,j,k) * ( A(k,row,col) * x(j,col) )
        if ( ( have_tensor = iA != kY && iA != jX ) ) { swap( kY , iA ); }
        break ;
      case 5 :
        // y(i,row) += tensor(i,j,k) * ( A(j,row,col) * x(k,col) )
        if (( have_tensor = iA != kY && iA != jX && kY != jX )) {
          swap(kY,iA); swap(iA,jX);
        }
        break ;
      }
    }

    // This thread is assigned:
    //   iBlockRow == block of the 'y' vector
    //   have_tensor, iA, jX, kY, tval == sparse tensor value
    // 
    for ( size_type iBlockRow = blockIdx.y ;
                    iBlockRow < m_A.graph.row_count() ;
                    iBlockRow += gridDim.y ) {

      VectorScalar y_sum = 0 ;

      const size_type iBlockEntryEnd = m_A.graph.row_entry_end(iBlockRow);
            size_type iBlockEntry    = m_A.graph.row_entry_begin(iBlockRow);

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        __syncthreads(); // Wait for previous query to complete

        if ( read_coeff ) {
          // Coalesced read of block for 'A' and 'x' into shared memory
          const size_type iBlockColumn = m_A.graph.column(iBlockEntry);

          shared_matrix_value( iCoeff ) = m_A.values( iCoeff , iBlockEntry );
          shared_vector_value( iCoeff ) = m_x( iCoeff , iBlockColumn );
        }

        __syncthreads(); // Wait for current load to complete

        if ( have_tensor ) {
          y_sum += tval * shared_matrix_value(iA) * shared_vector_value(jX);
        }
      }

      if ( have_tensor ) {
        cuda_internal_atomic_add( m_y( kY , iBlockRow ) , y_sum );
      }
    }
  }

};

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_CUDA_PRODUCTTENSOR_HPP */

