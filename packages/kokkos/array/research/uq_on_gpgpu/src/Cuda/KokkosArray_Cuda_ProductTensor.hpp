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

#ifndef KOKKOSARRAY_CUDA_PRODUCTTENSOR_HPP
#define KOKKOSARRAY_CUDA_PRODUCTTENSOR_HPP

#include <utility>
#include <sstream>
#include <stdexcept>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>

namespace KokkosArray {
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
    const size_type n = cuda_internal_maximum_warp_count() * CudaTraits::WarpSize ;

    if ( 3 * ( n / 3 ) < block.dimension() ) {
      throw std::runtime_error( std::string("KokkosArray::Impl::Multiply< SparseProductTensor<3,Scalar,Cuda> > ERROR: block too large") );
    }

    const dim3 d( std::min( block.entry_count() , n / 3 ) , 3 , 1 );

    return d ;
  }

  template< typename VectorValue >
  __host__
  static size_type shmem_size( const block_type & block )
  {
    const dim3 dim = thread_block( block );

    return 3 * block.dimension() * sizeof(VectorValue) +
           3 * dim.x * sizeof(int);
  }

  //----------------------------------------

  inline
  __device__
  static void swap( int & i , int & j )
  { const int tmp = i ; i = j ; j = tmp ; }

  //----------------------------------------

  template< typename MatrixValue , typename VectorValue >
  __device__
  static VectorValue apply( const block_type & block ,
                            const MatrixValue * const a ,
                            const VectorValue * const x )
  {
    const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;
    const size_type dimension = block.dimension();
    VectorValue * const shY = kokkos_impl_cuda_shared_memory<VectorValue>();
    VectorValue * const shX = (VectorValue *)( shY + dimension );
    VectorValue * const shA = (VectorValue *)( shX + dimension );

#if 0
    int * const shC = (int *)( shA + blockDim.x );
#endif

    // Initialize output, load shared input

    if ( tid < dimension ) {
      shY[ tid ] = 0 ;
      shX[ tid ] = x[ tid ];
      shA[ tid ] = a[ tid ];
    }

    __syncthreads(); // wait for all data reads

    // Loop through sparse tensor contributions
    // with batch reads.

    const size_type nEntryLoop = block.entry_count() + threadIdx.x ;

    for ( size_type iEntry = threadIdx.x ;
                    iEntry < nEntryLoop ; iEntry += blockDim.x ) {

#if 0

      shC[ tid ] = block.coord( iEntry , threadIdx.y );

      __syncthreads();

      if ( iEntry < block.entry_count() ) {
        int i = shC[ threadIdx.x ];
        int j = shC[ threadIdx.x + blockDim.x ];
        int k = shC[ threadIdx.x + blockDim.x * 2 ];

#else

      if ( iEntry < block.entry_count() ) {

        int i = block.coord( iEntry , 0 ); // coalesced
        int j = block.coord( iEntry , 1 ); // coalesced
        int k = block.coord( iEntry , 2 ); // coalesced
#endif

        if      ( 0 == threadIdx.y ) { }
        else if ( 1 == threadIdx.y && j != k ) { swap( k , j ); }
        else if ( 2 == threadIdx.y && i != k && i != j ) { swap( k , i ); }
        else { k = -1 ; }

        if ( 0 <= k ) {
          VectorValue tmp = shA[i] * shX[j] ;
          if ( i != j ) tmp += shA[j] * shX[i] ;
          cuda_internal_atomic_add( shY[k] , tmp * block.value( iEntry ) );
        }
      }

#if 0
      __syncthreads();
#endif
    }

    __syncthreads();

    return tid < dimension ? shY[ tid ] : 0 ;
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
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef SparseProductTensor< 3 , TensorScalar , device_type >    tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >           vector_type ;

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
      msg << "KokkosArray::Multiply< BlockCrsMatrix< SparseProductTensor<> ... > > "
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

    const unsigned row_count = A.graph.row_map.dimension(0) - 1 ;

    const size_type nGridY =
       std::min( row_count , Impl::CudaTraits::UpperBoundGridCount / nGridX );

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
  void operator()(void) const
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

    const unsigned row_count = m_A.graph.row_map.dimension(0) - 1 ;

    // This thread is assigned:
    //   iBlockRow == block of the 'y' vector
    //   have_tensor, iA, jX, kY, tval == sparse tensor value
    // 
    for ( size_type iBlockRow = blockIdx.y ;
                    iBlockRow < row_count ;
                    iBlockRow += gridDim.y ) {

      VectorScalar y_sum = 0 ;

      const size_type iBlockEntryEnd = m_A.graph.row_map[iBlockRow+1];
            size_type iBlockEntry    = m_A.graph.row_map[iBlockRow];

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        __syncthreads(); // Wait for previous query to complete

        if ( read_coeff ) {
          // Coalesced read of block for 'A' and 'x' into shared memory
          const size_type iBlockColumn = m_A.graph.entries(iBlockEntry);

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

#if 1

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< 3 , TensorScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensor< 3 , TensorScalar , device_type >       tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >           vector_type ;

  // m_A.graph.row_map.length() == gridDim.y
  // m_A.block.dimension() <= gridDim.x * blockDim.y

  class ProductTensorLoop {
  public:

    const matrix_type m_A ;
    const vector_type m_x ;
    const vector_type m_y ;

    ProductTensorLoop( const matrix_type & A ,
                       const vector_type & x ,
                       const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();

      VectorScalar * const sh_x = kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_A = sh_x + dim ;

      const size_type nid = CudaTraits::WarpSize * blockDim.y ;
      const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

      // This block is responsible for computing values of 'y'
      // in the range [ y_begin , y_end )
      //   y_begin = blockIdx.x * blockDim.y
      //   y_end   = y_begin + blockDim.y
      //
      // Value of 'y' for which this warp is responsible:
      //
      const size_type iyInner = threadIdx.y + blockDim.y * blockIdx.x ;

      // blockIdx.y == row in the deterministic (finite element) system
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.y + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ blockIdx.y ];

      // Product tensor entries which this warp will iterate:
      //
      const size_type iBeg = iyInner < dim ? m_A.block.entry_begin( iyInner ) : 0 ;
      const size_type iEnd = iyInner < dim ? m_A.block.entry_end(   iyInner ) : 0 ;

      VectorScalar y = 0 ;

      // Loop over columns in the discrete (finite element) system.

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );

        // Coalesced read of X and A into shared memory

        for ( size_type i = tid ; i < dim ; i += nid ) {
          sh_x[i] = m_x( i , iBlockColumn );
        }

        for ( size_type i = tid ; i < dim ; i += nid ) {
          sh_A[i] = m_A.values( i , iBlockEntry );
        }

        __syncthreads(); // wait for X and A to be read

        // Loop through sparse tensor contributions with coalesced reads.

        for ( size_type i = iBeg + threadIdx.x ;
                        i < iEnd ; i += CudaTraits::WarpSize ) {

          // Read 'CudaTraits::WarpSize' entries from the tensor
          const int j = m_A.block.coord( i , 0 ); // coalesced read
          const int k = m_A.block.coord( i , 1 ); // coalesced read

          VectorScalar tmp = sh_A[j] * sh_x[k] ;
          if ( j != k ) tmp += sh_A[k] * sh_x[j] ;

          y += tmp * m_A.block.value(i); // coalesced read
        }

        __syncthreads(); // wait for X and A to be used.
      }

      // Reduction of 'y' within 'CudaTraits::WarpSize'

      sh_x[ tid ] = y ;

      if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+16];
      if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 8];
      if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 4];
      if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 2];
      if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 1];

      if ( iyInner < dim && 0 == threadIdx.x ) {
        m_y( iyInner , blockIdx.y ) = sh_x[ tid ];
      }
    }
  };

  //------------------------------------
  // m_A.graph.row_map.length() == gridDim.y
  // m_A.block.dimension() <= gridDim.x * blockDim.y

  class ProductTensorOnce {
  public:

    const matrix_type m_A ;
    const vector_type m_x ;
    const vector_type m_y ;

    ProductTensorOnce( const matrix_type & A ,
                       const vector_type & x ,
                       const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      volatile VectorScalar * const sh =
        kokkos_impl_cuda_shared_memory<VectorScalar>();

      const size_type dim = m_A.block.dimension();
      const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;
      const size_type nid = CudaTraits::WarpSize * blockDim.y ;

      // Value of 'y' for which this thread is responsible
      const size_type iyInner = threadIdx.y + blockDim.y * blockIdx.x ;

      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.y + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ blockIdx.y ];

      int j = -1 ;
      int k = -1 ;
      VectorScalar v = 0 ;

      if ( iyInner < dim ) {
        const size_type iP = threadIdx.x + m_A.block.entry_begin( iyInner );
        if ( iP < m_A.block.entry_end( iyInner ) ) {
          j = m_A.block.coord( iP , 0 );
          k = m_A.block.coord( iP , 1 );
          v = m_A.block.value( iP );
        }
      }

      VectorScalar y = 0 ;

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        // Coalesced read of X and A into shared memory

        for ( size_type i = tid ; i < dim ; i += nid ) {
          sh[ i ] = m_x( i , m_A.graph.entries(iBlockEntry) );
        }

        for ( size_type i = tid ; i < dim ; i += nid ) {
          sh[ i + dim ] = m_A.values( i , iBlockEntry );
        }

        __syncthreads(); // wait for X and A to be read

        if ( 0 <= j ) {
          VectorScalar tmp = sh[j+dim] * sh[k] ;
          if ( j != k ) tmp += sh[k+dim] * sh[j] ;

          y += tmp * v ;
        }

        __syncthreads(); // wait for X and A to be used.
      }

      // Reduction of 'y' within 'CudaTraits::WarpSize'

      sh[ tid ] = y ;

      if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh[tid] += sh[tid + 16];
      if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh[tid] += sh[tid +  8];
      if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh[tid] += sh[tid +  4];
      if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh[tid] += sh[tid +  2];
      if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh[tid] += sh[tid +  1];

      if ( iyInner < dim && 0 == threadIdx.x ) {
        m_y( iyInner , blockIdx.y ) = sh[ tid ];
      }
    }
  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    // Have at least nPool warps working on each tensor-block

    const size_type row_count = A.graph.row_map.dimension(0) - 1 ;

    const size_type nPool = 8 ;

    const size_type maxWarp =
      std::min( (A.block.dimension()+nPool-1) / nPool ,
                cuda_internal_maximum_warp_count() );

    size_type nWarp = 2 ; // need at least two warps

    for ( ; ( nWarp << 1 ) <= maxWarp ; nWarp <<= 1 );

    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );

    // Fill out block dimension with X
    // Fill out graph block-row count with Y
    const dim3 dGrid( ( A.block.dimension() + dBlock.y - 1 ) / dBlock.y ,
                      row_count , 1 );

    // Maximum required to:
    //   Reduce one value of Y per thread
    //   Read block of A and X into shared memory
    const size_type shmem =
      sizeof(VectorScalar) * std::max( dBlock.x * dBlock.y ,
                                       2 * A.block.dimension() );

    if ( 0 && A.block.entry_maximum() <= dBlock.x ) {
      // Disable this option for now ...
      // Product tensor entries can be read and held by a thread
      Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
        ( ProductTensorOnce( A , x , y ) );
    }
    else {
      // Must loop through product tensor entries.
      Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
        ( ProductTensorLoop( A , x , y ) );
    }
  }
};

#elif 0

// Maximum shared memory use variant:

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< 3 , TensorScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensor< 3 , TensorScalar , device_type >       tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >           vector_type ;

  // m_A.graph.row_map.length() == gridDim.y
  // m_A.block.dimension() <= gridDim.x * blockDim.y

  class ProductTensorLoop {
  public:

    const matrix_type m_A ;
    const vector_type m_x ;
    const vector_type m_y ;

    ProductTensorLoop( const matrix_type & A ,
                       const vector_type & x ,
                       const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type max = m_A.block.entry_maximum();

      // Shared memory:
      //   sh_x[ tensor_dimension ]
      //   sh_A[ tensor_dimension ]
      //
      //   sh_block_T[ tensor_entry_max * nWarp ]
      //   sh_block_j[ tensor_entry_max * nWarp ]
      //   sh_block_k[ tensor_entry_max * nWarp ]
      //

      VectorScalar * const sh_x = kokkos_impl_cuda_shared_memory<VectorScalar>();
      MatrixScalar * const sh_A = (MatrixScalar *)( sh_x + dim );
      TensorScalar * const sh_block_T = (TensorScalar *)( sh_A + dim );
      size_type    * const sh_block_j = (size_type *)( sh_block_T + max * blockDim.y );
      size_type    * const sh_block_k = (size_type *)( sh_block_j + max * blockDim.y );

      // This warp's shared memory for tensor data:

      TensorScalar * const sh_warp_T = sh_block_T + max * threadIdx.y ;
      size_type    * const sh_warp_j = sh_block_j + max * threadIdx.y ;
      size_type    * const sh_warp_k = sh_block_k + max * threadIdx.y ;

      // Block size and thread id within the entire block:

      const size_type nid = CudaTraits::WarpSize * blockDim.y ;
      const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

      // This block is responsible for computing values of 'y'
      // in the range [ y_begin , y_end )
      //   y_begin = blockIdx.x * blockDim.y
      //   y_end   = y_begin + blockDim.y
      //
      // Value of 'y' for which this warp is responsible:
      //
      const size_type iyInner = threadIdx.y + blockDim.y * blockIdx.x ;

      // Product tensor entries which this warp will iterate:
      //
      const size_type iBeg = iyInner < dim ? m_A.block.entry_begin( iyInner ) : 0 ;
      const size_type iEnd = iyInner < dim ? m_A.block.entry_end(   iyInner ) : 0 ;
      const size_type iCount = iEnd - iBeg ;

      // Each warp reads its own tensor data into shared memory:

      for ( size_type i = threadIdx.x ;
                      i < iCount ; i += CudaTraits::WarpSize ) {
        sh_warp_j[i] = m_A.block.coord( i + iBeg , 0 ); // coalesced read
        sh_warp_k[i] = m_A.block.coord( i + iBeg , 1 ); // coalesced read
        sh_warp_T[i] = m_A.block.value( i + iBeg );     // coalesced read
      }

      // blockIdx.y == row in the deterministic (finite element) system
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.y + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ blockIdx.y ];

      VectorScalar y = 0 ;

      // Loop over columns in the discrete (finite element) system.

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );
        const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

        // Whole block reads X and A into shared memory:

        for ( size_type i = tid ; i < dim ; i += nid ) {
          sh_x[i] = x[i] ; // m_x(        i , iBlockColumn );
          sh_A[i] = A[i] ; // m_A.values( i , iBlockEntry );
        }

        __syncthreads(); // wait for X and A to be read

        // Loop through sparse tensor contributions:

        for ( size_type i = threadIdx.x ;
                        i < iCount ; i += CudaTraits::WarpSize ) {

          const int j = sh_warp_j[i] ;
          const int k = sh_warp_k[i] ;

          VectorScalar tmp = sh_A[j] * sh_x[k] ;

          if ( j != k ) { tmp += sh_A[k] * sh_x[j] ; }

          y += tmp * sh_warp_T[i] ;
        }

        __syncthreads(); // wait for X and A to be used.
      }

      // Reduction of 'y' within 'CudaTraits::WarpSize'

      sh_x[ tid ] = y ;

      if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+16];
      if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 8];
      if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 4];
      if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 2];
      if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_x[tid] += sh_x[tid+ 1];

      if ( iyInner < dim && 0 == threadIdx.x ) {
        m_y( iyInner , blockIdx.y ) = sh_x[ tid ];
      }
    }
  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count        = A.graph.row_map.dimension(0) - 1 ;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tensor_entry_max = A.block.entry_maximum();

    // Each warp computes a single value of 'y'
    //
    // Shared memory:
    //   sh_x[ tensor_dimension ]
    //   sh_A[ tensor_dimension ]
    //
    //   sh_block_T[ tensor_entry_max * nWarp ]
    //   sh_block_j[ tensor_entry_max * nWarp ]
    //   sh_block_k[ tensor_entry_max * nWarp ]
    //
    const int nsh_base = sizeof(VectorScalar) * tensor_dimension +
                         sizeof(MatrixScalar) * tensor_dimension ;

    const int nsh_warp = sizeof(TensorScalar) * tensor_entry_max +
                         sizeof(size_type)    * tensor_entry_max * 2 ;

    enum { ShCapacity = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity / 2 };

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

    // Must have at least one VectorScalar per thread:

    const size_type shmem = std::max( (size_type)( nsh_base + nsh_warp * nWarp ),
                                      (size_type)( sizeof(VectorScalar) * CudaTraits::WarpSize * nWarp ) );

    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );

    // Fill out block dimension with X
    // Fill out graph block-row count with Y
    const dim3 dGrid( ( tensor_dimension + nWarp - 1 ) / nWarp ,
                      row_count , 1 );

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem << ")" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_entry_max(" << tensor_entry_max << ")" << std::endl
              ;
#endif

    Impl::CudaParallelLaunch< ProductTensorLoop >( ProductTensorLoop( A , x , y ), dGrid , dBlock , shmem );
  }
};

#else

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< 3 , TensorScalar , Cuda > ,
                  MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensor< 3 , TensorScalar , device_type >       tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >           vector_type ;

  // m_A.graph.row_map.length() == gridDim.x

  class ProductTensorLoop {
  public:

    const matrix_type m_A ;
    const vector_type m_x ;
    const vector_type m_y ;

    ProductTensorLoop( const matrix_type & A ,
                       const vector_type & x ,
                       const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type tensor_dim = m_A.block.dimension();
      const size_type tensor_dim_align = KokkosArray::Impl::CudaTraits::warp_align( tensor_dim );

      // Shared memory:
      //   sh_x[ tensor_dim_align ]
      //   sh_A[ tensor_dim_align ]
      //   sh_y[ CudaTraits::WarpSize * blockDim.y ]

      VectorScalar * const sh_x = kokkos_impl_cuda_shared_memory<VectorScalar>();
      MatrixScalar * const sh_A = (MatrixScalar *)( sh_x + tensor_dim_align );
      volatile VectorScalar * const sh_y = (VectorScalar *)( sh_A + tensor_dim_align );
      unsigned int * const sh_work = (unsigned int *)( sh_y + CudaTraits::WarpSize * blockDim.y );

      // Block size and thread id within the entire block:

      const size_type nid = CudaTraits::WarpSize * blockDim.y ;
      const size_type tid = threadIdx.x + CudaTraits::WarpSize * threadIdx.y ;

      // Loop over columns in the discrete (finite element) system.

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
            size_type iBlockEntry    = m_A.graph.row_map[ blockIdx.x ];

      for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {

        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );

        const VectorScalar * const x = & m_x(        0 , iBlockColumn );
        const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

        // Whole block reads X and A into shared memory:

        for ( size_type i = tid ; i < tensor_dim ; i += nid ) {
          sh_x[i] = x[i] ; // m_x(        i , iBlockColumn );
          sh_A[i] = A[i] ; // m_A.values( i , iBlockEntry );
        }

        if ( 0 == tid ) { *sh_work = 0 ; } // For work stealing strategy

        __syncthreads(); // wait for X and A to be read and sh_work to be set

        // Each warp computes a row of the stochastic block
        // for coalesced reads and no need for explicit synchronization.

#if 1
        // Simple striding
        for ( size_type iyInner = threadIdx.y ; iyInner < tensor_dim ; iyInner += blockDim.y ) {

#else
        // Try work stealing since rows are irregular length
        for (;;) {

          if ( 0 == threadIdx.x ) { sh_work[ threadIdx.y + 1 ] = atomicInc( sh_work , ~0u ); }

          const unsigned int iyInner = sh_work[ threadIdx.y + 1 ];

          if ( tensor_dim <= iyInner ) break ;
#endif

          VectorScalar y = 0 ;

          // Product tensor entries which this warp will iterate:
          //
          const size_type iBeg = m_A.block.entry_begin( iyInner );
          const size_type iEnd = m_A.block.entry_end(   iyInner );


          // Loop through sparse tensor contributions:

          for ( size_type i = threadIdx.x + iBeg ;
                          i < iEnd ; i += CudaTraits::WarpSize ) {

            const int j = m_A.block.coord( i , 0 ); // coalesced read
            const int k = m_A.block.coord( i , 1 ); // coalesced read

            y += m_A.block.value( i ) * ( sh_A[j] * sh_x[k] + ( j != k ? sh_A[k] * sh_x[j] : 0 ) );
          }

          // Reduction of 'y' within the warp

          sh_y[ tid ] = y ;

          if ( threadIdx.x + 16 < CudaTraits::WarpSize ) sh_y[tid] += sh_y[tid+16];
          if ( threadIdx.x +  8 < CudaTraits::WarpSize ) sh_y[tid] += sh_y[tid+ 8];
          if ( threadIdx.x +  4 < CudaTraits::WarpSize ) sh_y[tid] += sh_y[tid+ 4];
          if ( threadIdx.x +  2 < CudaTraits::WarpSize ) sh_y[tid] += sh_y[tid+ 2];
          if ( threadIdx.x +  1 < CudaTraits::WarpSize ) sh_y[tid] += sh_y[tid+ 1];

          if ( 0 == threadIdx.x ) {
            m_y( iyInner , blockIdx.x ) = sh_y[ tid ];
          }
        }

        __syncthreads(); // wait for X and A to be used before updating
      }
    }
  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count        = A.graph.row_map.dimension(0) - 1 ;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tensor_dim_align = KokkosArray::Impl::CudaTraits::warp_align( tensor_dimension );

    // Shared memory:
    //   sh_x[ tensor_dim_align ]
    //   sh_A[ tensor_dim_align ]
    //   sh_y[ CudaTraits::WarpSize * blockDim.y ]
    //   sh_work[ 1 + blockDim.y ]
    //
    // tensor_dim <= m_sh_stride && m_sh_stride % WarpSize == 0

    const int nsh_base = sizeof(VectorScalar) * tensor_dim_align +
                         sizeof(MatrixScalar) * tensor_dim_align +
                         sizeof(unsigned int) ;

    const int nsh_warp = sizeof(VectorScalar) * CudaTraits::WarpSize +
                         sizeof(unsigned int) ;

    enum { ShCapacity = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity / 2 };

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

    // Must have at least one VectorScalar per thread:

    const size_type shmem = nsh_base + nsh_warp * nWarp + sizeof(unsigned int);

    const dim3 dBlock( CudaTraits::WarpSize , nWarp , 1 );

    const dim3 dGrid( row_count , 1 , 1 );

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem << ")" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_dim_align(" << tensor_dim_align << ")" << std::endl
              ;
#endif

    cuda_parallel_launch_local_memory< ProductTensorLoop ><<< dGrid , dBlock , shmem >>>( ProductTensorLoop( A , x , y ) );
    // Impl::CudaParallelLaunch< ProductTensorLoop >( ProductTensorLoop( A , x , y ), dGrid , dBlock , shmem );
  }
};

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_CUDA_PRODUCTTENSOR_HPP */

