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
//----------------------------------------------------------------------------

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


  static void execute( const matrix_type & A ,
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

    if ( nWarp  > Impl::cuda_internal_maximum_warp_count() ||
         nGridX > Impl::CudaTraits::UpperBoundGridCount ||
         shared_size * sizeof(size_type) >
            Impl::cuda_internal_maximum_shared_words() ) {

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
      if ( shared_size * sizeof(size_type) >
             Impl::cuda_internal_maximum_shared_words() ) {
        msg << " : shared_memory " << shared_size << " > " 
            <<  Impl::cuda_internal_maximum_shared_words() * sizeof(size_type);
      }
      throw std::runtime_error( msg.str() );
    }

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

      const size_type iBlockEntryEnd = m_A.graph.row_range_end(iBlockRow);
            size_type iBlockEntry    = m_A.graph.row_range_begin(iBlockRow);

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
        atomicAdd( & m_y( kY , iBlockRow ) , y_sum );
      }
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_CUDA_PRODUCTTENSOR_HPP */

