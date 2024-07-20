// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_CRSMATRIX_UQ_PCE_CUDA_HPP
#define KOKKOS_CRSMATRIX_UQ_PCE_CUDA_HPP

#if defined( __CUDACC__)

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

//MD 08/2017 Note: I commented out below, as this file is totally
//removed from KokkosKernels. It does not look like the
//included file is used anywhere in the file.
//#include "Kokkos_MV.hpp" // for some utilities

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsProductTensor.hpp"

#include "Kokkos_Core.hpp"

#include "Stokhos_Cuda_DeviceProp.hpp"
//#include "Stokhos_Cuda_WarpShuffle.hpp"

#include "Teuchos_TestForException.hpp"

//#include "cuda_profiler_api.h"

#ifdef __CUDA_ARCH__
#  if (__CUDA_ARCH__ >= 300)
#    define HAVE_CUDA_SHUFFLE 1
#  else
#    define HAVE_CUDA_SHUFFLE 0
#  endif
#else
#  define HAVE_CUDA_SHUFFLE 0
#endif

namespace Stokhos {

//----------------------------------------------------------------------------
// Specialization of Kokkos CrsMatrix math functions
//----------------------------------------------------------------------------

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                   MatrixOrdinal,
                                   Kokkos::Cuda,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< const Sacado::UQ::PCE<InputStorage>*,
                              InputP... >,
                Kokkos::View< Sacado::UQ::PCE<OutputStorage>*,
                              OutputP... >
                >
{
public:
  typedef Sacado::UQ::PCE<MatrixStorage> MatrixValue;
  typedef Sacado::UQ::PCE<InputStorage> InputVectorValue;
  typedef Sacado::UQ::PCE<OutputStorage> OutputVectorValue;

  typedef Kokkos::Cuda MatrixDevice;
  typedef MatrixDevice execution_space;
  typedef execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                             MatrixOrdinal,
                             MatrixDevice,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename Kokkos::CijkType<matrix_values_type>::type tensor_type;
  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;

private:

  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename matrix_values_type::array_type matrix_array_type;
  typedef typename input_vector_type::array_type input_array_type;
  typedef typename output_vector_type::array_type output_array_type;

  typedef typename MatrixValue::value_type matrix_scalar;
  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;
  typedef typename tensor_type::value_type tensor_scalar;

  const matrix_array_type   m_A_values ;
  const matrix_graph_type   m_A_graph ;
  const input_array_type    m_x ;
  const output_array_type   m_y ;
  const tensor_type         m_tensor ;
  const input_scalar        m_a ;
  const output_scalar       m_b ;
  const size_type           BlockSize;

  Multiply( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y ,
            const input_scalar & a ,
            const output_scalar & b ,
            const size_type block_size )
  : m_A_values( A.values )
  , m_A_graph( A.graph )
  , m_x( x )
  , m_y( y )
  , m_tensor( Kokkos::cijk(A.values) )
  , m_a( a )
  , m_b( b )
  , BlockSize(block_size)
  {}

public:

  __device__ void operator()(void) const
  {
    // Number of bases in the stochastic system:
    const size_type dim = m_tensor.dimension();

    // Get shared memory for loading x, A, and y
    volatile input_scalar * const sh_x =
      kokkos_impl_cuda_shared_memory<input_scalar>();
    volatile matrix_scalar * const sh_A = sh_x + BlockSize*dim;
    volatile output_scalar * const sh_y = sh_A + BlockSize*dim;
#if !HAVE_CUDA_SHUFFLE
    volatile output_scalar * const sh_t = sh_y + dim;
#endif

    const size_type nid = blockDim.x * blockDim.y;
    const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

    // Zero y
    for ( size_type i = tid; i < dim; i += nid ) {
      sh_y[i] = 0.0;
    }

    // Loop over columns in the discrete (finite element) system.
    // blockIdx.x == row in the deterministic (finite element) system
    const size_type iBlockEntryBeg = m_A_graph.row_map[ blockIdx.x ];
    const size_type iBlockEntryEnd = m_A_graph.row_map[ blockIdx.x + 1 ];
    for (size_type iBlockEntry=iBlockEntryBeg; iBlockEntry<iBlockEntryEnd;
         iBlockEntry += BlockSize) {
      const size_type block_size =
        (iBlockEntryEnd-iBlockEntry < BlockSize) ?
        iBlockEntryEnd-iBlockEntry : BlockSize;

      // Wait for X and A to be used in the previous iteration
      // before reading new values.
      __syncthreads();

      // Coalesced read blocks of X and A into shared memory
      for ( size_type col = 0; col < block_size; ++col ) {

        const size_type iBlockColumn = m_A_graph.entries( iBlockEntry + col );
        const input_scalar * const x =  & m_x(        0, iBlockColumn );
        const matrix_scalar * const A = & m_A_values( iBlockEntry + col, 0 );

        // Coalesced read by the whole block from global memory:
        for ( size_type i = tid; i < dim; i += nid ) {
          sh_x[col + i * BlockSize] = x[i]; // m_x(        i, iBlockColumn );
          sh_A[col + i * BlockSize] = A[i]; // m_A_values( iBlockEntry , i );
        }

      }

      __syncthreads(); // wait for X and A to be read before using them

      // This cuda block is responsible for computing all values of 'y'
      for ( size_type i = threadIdx.y; i < dim; i += blockDim.y ) {
        output_scalar y = 0;

        // Product tensor entries which this warp will iterate:
        const size_type lBeg = m_tensor.entry_begin( i );
        const size_type lEnd = m_tensor.entry_end(   i );

        // Loop through sparse tensor contributions with coalesced reads.
        for ( size_type l = lBeg+threadIdx.x; l < lEnd; l += blockDim.x ) {

          // Read 'blockDim.x' entries from the tensor (coalesced)
          const size_type kj    = m_tensor.coord( l );
          const tensor_scalar v = m_tensor.value( l );
          const size_type j     = ( kj & 0x0ffff ) * BlockSize ;
          const size_type k     = ( kj >> 16     ) * BlockSize ;

          for ( size_type col = 0; col < block_size; ++col ) {
            y += v * ( sh_A[col+j] * sh_x[col+k] +
                       sh_A[col+k] * sh_x[col+j] );
          }

        }

        // Reduction of 'y' within 'blockDim.x'
#if HAVE_CUDA_SHUFFLE
        if (blockDim.x >= 2)  y += Kokkos::shfl_down(y, 1,  blockDim.x);
        if (blockDim.x >= 4)  y += Kokkos::shfl_down(y, 2,  blockDim.x);
        if (blockDim.x >= 8)  y += Kokkos::shfl_down(y, 4,  blockDim.x);
        if (blockDim.x >= 16) y += Kokkos::shfl_down(y, 8,  blockDim.x);
        if (blockDim.x >= 32) y += Kokkos::shfl_down(y, 16, blockDim.x);
        if ( threadIdx.x == 0 ) sh_y[i] += y;
#else
        sh_t[ tid ] = y;
        if (threadIdx.x+16 < blockDim.x) sh_t[tid] += sh_t[tid+16];
        if (threadIdx.x+ 8 < blockDim.x) sh_t[tid] += sh_t[tid+ 8];
        if (threadIdx.x+ 4 < blockDim.x) sh_t[tid] += sh_t[tid+ 4];
        if (threadIdx.x+ 2 < blockDim.x) sh_t[tid] += sh_t[tid+ 2];
        if (threadIdx.x+ 1 < blockDim.x) sh_t[tid] += sh_t[tid+ 1];
        if (threadIdx.x == 0) sh_y[i] += sh_t[tid];
#endif

      }

    }

    // Wait for all contributions of y to be completed
    __syncthreads();

    // Store result back in global memory
    if ( m_b == output_scalar(0) )
      for ( size_type i = tid; i < dim; i += nid )
        m_y( i, blockIdx.x ) = m_a * sh_y[ i ];
    else
      for ( size_type i = tid; i < dim; i += nid )
        m_y( i, blockIdx.x ) = m_a * sh_y[ i ] + m_b * m_y( i, blockIdx.x );
  }

  struct TensorReadEntry {
    size_type block_size, shmem, num_blocks, num_warp;
    double reads;
  };

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const tensor_type tensor = Kokkos::cijk(A.values);
    const size_type row_count = A.graph.row_map.extent(0) - 1;
    const size_type tensor_dimension = tensor.dimension();
    const size_type tensor_align = tensor_dimension;
    const size_type avg_tensor_entries_per_row = tensor.avg_entries_per_row();

    // Should compute this from FEM graph
    const size_type fem_nnz_per_row = 27;

    // Get device properties we need for whatever device is currently selected
    DeviceProp device_prop;
    const size_type shcap = device_prop.shared_memory_capacity;
    const size_type sh_granularity = device_prop.shared_memory_granularity;
    const size_type max_shmem_per_block = device_prop.max_shmem_per_block;
    const size_type max_blocks_per_sm = device_prop.max_blocks_per_sm;
    const size_type warp_size = device_prop.warp_size;
    const size_type warp_granularity = device_prop.warp_granularity;
    const size_type max_warps_per_block =
      std::min(device_prop.max_threads_per_block / warp_size,
               device_prop.max_warps_per_sm);
    const size_type min_warps_per_block = 1;
    const size_type max_regs_per_sm = device_prop.max_regs_per_sm;
    const size_type max_regs_per_block = device_prop.max_regs_per_block;
    const size_type reg_bank_size = device_prop.reg_bank_size;

    // Compute number of warps we can fit on each SM based on register limits
    // Use Cuda introspection to determine number of registers per thread
    //const size_type regs_per_thread = 46;
    const size_type regs_per_thread =
      device_prop.get_kernel_registers(
        Kokkos::Impl::cuda_parallel_launch_local_memory<Multiply>);
    const size_type regs_per_warp =
      (warp_size*regs_per_thread + reg_bank_size-1) & ~(reg_bank_size-1);
    const size_type warps_per_sm =
      (max_regs_per_sm/regs_per_warp) & ~(warp_granularity-1);
    const size_type warps_per_block =
      (max_regs_per_block/regs_per_warp) & ~(warp_granularity-1);

    // Compute number of threads per stochastic row based on number of
    // nonzero entries per row.
    // For double, 16 threads/row is still coalesced, but not for float.
    // We should reorder the tensor values for a given vector width to
    // maintain coalesced reads.  This would help smaller problems too by
    // allowing fewer threads/row.
    const size_type threads_per_row =
      avg_tensor_entries_per_row >= 88 ? 32 : 16;
    const size_type rows_per_warp = warp_size / threads_per_row;

    const size_type in_vec_scalar_size = sizeof(input_scalar);
    const size_type out_vec_scalar_size = sizeof(output_scalar);
    const size_type mat_scalar_size = sizeof(matrix_scalar);

#define USE_FIXED_BLOCKSIZE 0

#if USE_FIXED_BLOCKSIZE

    const size_type num_blocks = 3;
    size_type nw = warps_per_sm / num_blocks;
    while (nw > 1 && num_blocks*nw % warp_granularity) --nw;
    const size_type num_warp = nw;
    const size_type sh_per_block = shcap / num_blocks;
    const size_type sr =
      device_prop.has_shuffle ? 0 : in_vec_scalar_size*warp_size*num_warp;
    size_type bs = ((sh_per_block - sr) / tensor_align - out_vec_scalar_size) /
      (in_vec_scalar_size+mat_scalar_size);
    if (bs % 2 == 0) --bs;
    const size_type block_size_max = 31;
    const size_type block_size = std::min(bs, block_size_max);
    //const size_type block_size = 7;
    const size_type shmem =
      ( ((in_vec_scalar_size+mat_scalar_size)*block_size+out_vec_scalar_size)*tensor_align + sr + sh_granularity-1 ) & ~(sh_granularity-1);

#else

    // We want to maximize the number of blocks per SM (to maximize throughput)
    // as well as the block_size (to minimize tensor reads), subject to
    // shared memory and register constraints.  Here we try to do this computing
    // the number of tensor reads per block per thread for each choice of
    // block_size, and then choose the configuration with the smallest value.
    // This isn't perfect, but seems to generally work OK.  It could be
    // improved with a better model of:
    //   * Number of blocks versus warps per block (to minimize synchronization)
    //   * Thread efficiency for small numbers of rows per thread
    const size_type block_size_min = 3;
    const size_type half_nnz_per_row = fem_nnz_per_row / 2 + 1;
    const size_type block_size_max =
      half_nnz_per_row % 2 ? half_nnz_per_row + 1 : half_nnz_per_row;
    Teuchos::Array<TensorReadEntry> reads_per_thread;
    for (size_type bs = block_size_min; bs<=block_size_max; bs+=2) {
      // We don't know the number of warps yet, so we just have to bound
      // sr by the maximum number possible (which is all warps in 1 block)
      const size_type sr =
        device_prop.has_shuffle ? 0 : in_vec_scalar_size*warp_size*warps_per_block;
      size_type shmem =
        ((in_vec_scalar_size+mat_scalar_size)*bs+out_vec_scalar_size)*tensor_align+sr;
      shmem = (shmem + sh_granularity-1) & ~(sh_granularity-1);
      if (shmem <= max_shmem_per_block) {
        size_type num_blocks = std::min(shcap / shmem, max_blocks_per_sm);
        size_type tensor_reads = (fem_nnz_per_row+bs-1) / bs;
        size_type num_warp =
          std::min(std::max(std::min(warps_per_sm/num_blocks, warps_per_block),
                            min_warps_per_block),
                   max_warps_per_block);
        while (num_warp > 1 && num_blocks*num_warp % warp_granularity)
          --num_warp;
        TensorReadEntry entry;
        entry.block_size = bs;
        entry.shmem = shmem;
        entry.num_blocks = num_blocks;
        entry.num_warp = num_warp;

        // Prefer at least 3 blocks
        size_type factor = std::min(num_blocks,3u);
        entry.reads = (static_cast<double>(tensor_reads) /
                       static_cast<double>(factor*num_blocks*num_warp));
        reads_per_thread.push_back(entry);
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      reads_per_thread.size() == 0, std::logic_error,
      "Stochastic problem dimension is too large to fit in shared memory");
    size_type idx = 0;
    double min_reads = reads_per_thread[0].reads;
    for (int i=1; i<reads_per_thread.size(); ++i) {
      if (reads_per_thread[i].reads < min_reads) {
        idx = i;
        min_reads = reads_per_thread[i].reads;
      }
    }

    const size_type block_size = reads_per_thread[idx].block_size;
    const size_type shmem = reads_per_thread[idx].shmem;
    const size_type num_blocks = reads_per_thread[idx].num_blocks;
    const size_type num_warp = reads_per_thread[idx].num_warp;

#endif

    // Setup thread blocks and grid
    const dim3 dBlock( threads_per_row , rows_per_warp*num_warp , 1 );
    const dim3 dGrid( row_count, 1, 1 );

#if 0
    std::cout << "block_size = " << block_size
              << " tensor reads = " << (fem_nnz_per_row+block_size-1)/block_size
              << " regs_per_thread = " << regs_per_thread
              << " num blocks = " << num_blocks
              << " num warps = " << num_warp
              << " num rows = " << tensor_dimension
              << " rows/warp = " << tensor_dimension / (num_warp*rows_per_warp)
              << " avg entries/row = " <<  avg_tensor_entries_per_row
              << std::endl;
#endif

    // Finally launch our kernel
    //cudaProfilerStart();
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< dGrid, dBlock, shmem >>>
      ( Multiply( A, x, y, a, b, block_size ) );
    //cudaProfilerStop();
  }
};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>**,...>,
//   x and y are rank 2
//
// Note:  Unlike the rank-1 version, this version has not been
// optimized, and doesn't even include the block-column implementation
template <typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                   MatrixOrdinal,
                                   Kokkos::Cuda,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< const Sacado::UQ::PCE<InputStorage>**,
                              InputP... >,
                Kokkos::View< Sacado::UQ::PCE<OutputStorage>**,
                              OutputP... >
                >
{
public:
  typedef Sacado::UQ::PCE<MatrixStorage> MatrixValue;
  typedef Sacado::UQ::PCE<InputStorage> InputVectorValue;
  typedef Sacado::UQ::PCE<OutputStorage> OutputVectorValue;

  typedef Kokkos::Cuda MatrixDevice;
  typedef MatrixDevice execution_space;
  typedef execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                             MatrixOrdinal,
                             MatrixDevice,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;
  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

public:

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    typedef Kokkos::View< const InputVectorValue*, InputP... > input_vector_type_1D;
    typedef Kokkos::View< OutputVectorValue*, OutputP... > output_vector_type_1D;
    typedef Multiply< matrix_type, input_vector_type_1D,
      output_vector_type_1D > multiply_type_1D;

    const size_type num_col = y.extent(1);
    for (size_type col=0; col<num_col; ++col)

      multiply_type_1D::apply(
        A,
        Kokkos::subview( x, Kokkos::ALL(), col),
        Kokkos::subview(y, Kokkos::ALL(), col),
        a, b );
  }
};

template <typename Kernel>
__global__ void
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1024,2)
#endif
MeanFullOccupancyKernelLaunch(Kernel kernel) {
  kernel();
}

// Kernel implementing y = A * x where PCE size of A is 1
//   A == Kokkos::CrsMatrix< Sacado::UQ::PCE<...>,...>, with A.values.sacado_size() == 1
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class MeanMultiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                       MatrixOrdinal,
                                       Kokkos::Cuda,
                                       MatrixMemory,
                                       MatrixSize >,
                    Kokkos::View< const Sacado::UQ::PCE<InputStorage>*,
                                  InputP... >,
                    Kokkos::View< Sacado::UQ::PCE<OutputStorage>*,
                                  OutputP... >,
                    >
{
public:
  typedef Sacado::UQ::PCE<MatrixStorage> MatrixValue;
  typedef Sacado::UQ::PCE<InputStorage> InputVectorValue;
  typedef Sacado::UQ::PCE<OutputStorage> OutputVectorValue;

  typedef Kokkos::Cuda MatrixDevice;
  typedef MatrixDevice execution_space;
  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                             MatrixOrdinal,
                             MatrixDevice,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename MatrixValue::ordinal_type size_type;
  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;

  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename MatrixValue::value_type matrix_scalar;
  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  template <int BlockSize>
  struct Kernel {
    typedef MatrixDevice execution_space;
    typedef typename Kokkos::FlatArrayType<matrix_values_type>::type matrix_array_type;
    typedef typename input_vector_type::array_type input_array_type;
    typedef typename output_vector_type::array_type output_array_type;

    const matrix_array_type   m_A_values ;
    const matrix_graph_type   m_A_graph ;
    const output_array_type   v_y ;
    const input_array_type    v_x ;
    const input_scalar        m_a ;
    const output_scalar       m_b ;
    const size_type           m_row_count;
    const size_type           dim ;

    Kernel( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y ,
            const input_scalar & a ,
            const output_scalar & b )
      : m_A_values( A.values )
      , m_A_graph( A.graph )
      , v_y( y )
      , v_x( x )
      , m_a( a )
      , m_b( b )
      , m_row_count( A.graph.row_map.extent(0)-1 )
      , dim( dimension_scalar(x) )
      {}

    __device__ void operator()(void) const
    {
      // Store matrix values and column indices in shared memory
      // to reduce global memory reads
      volatile matrix_scalar * const sh_A =
        kokkos_impl_cuda_shared_memory<matrix_scalar>();
      volatile size_type * const sh_col =
        reinterpret_cast<volatile size_type*>(sh_A + BlockSize*blockDim.y);

      // Check for valid row
      const size_type iBlockRow = blockDim.y*blockIdx.x + threadIdx.y;
      if (iBlockRow < m_row_count) {

        const size_type iEntryBegin = m_A_graph.row_map[ iBlockRow ];
        const size_type iEntryEnd   = m_A_graph.row_map[ iBlockRow + 1 ];

        // Initialize result
        if (m_b == output_scalar(0))
          for ( size_type pce = threadIdx.x; pce < dim ; pce+=blockDim.x )
            v_y(pce, iBlockRow) = 0.0;
        else
          for ( size_type pce = threadIdx.x; pce < dim ; pce+=blockDim.x )
            v_y(pce, iBlockRow) = m_b*v_y(pce, iBlockRow);

        // Loop over columns in chunks of size BlockSize
        for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
             col_block+=BlockSize) {
          const size_type num_col = col_block+BlockSize <= iEntryEnd ?
            BlockSize : iEntryEnd-col_block;

          // Read BlockSize entries column indices at a time to maintain
          // coalesced accesses
          for (size_type col=threadIdx.x; col<num_col; col+=blockDim.x) {
            sh_col[col*blockDim.y+threadIdx.y] =
              m_A_graph.entries(col_block+col);
            sh_A[col*blockDim.y+threadIdx.y] =
              m_A_values(col_block+col);
          }
          if (blockDim.x > Kokkos::Impl::CudaTraits::WarpSize)
            __syncthreads();

          // Do portion mat-vec for each PCE coefficient and for columns
          // within this block
          for ( size_type pce = threadIdx.x; pce < dim ; pce+=blockDim.x ) {
            output_scalar s = 0.0;
            for ( size_type col = 0; col < num_col; ++col ) {
              const size_type iCol = sh_col[col*blockDim.y+threadIdx.y];
              const matrix_scalar aA = m_a*sh_A[col*blockDim.y+threadIdx.y];
              s += aA*v_x(pce, iCol);
            }
            v_y(pce, iBlockRow) += s;
          }
        }
      }
    }
  };

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const size_t row_count = A.graph.row_map.extent(0) - 1;
    const size_type dim = dimension_scalar(x);

    // Compute number of threads for PCE coefficients and number of
    // matrix rows processed by each CUDA block.  A total of 256 threads
    // gives best occupancy
    size_type threads_per_row;
    size_type rows_per_block;
    if (dim >= 32) {
      threads_per_row = 32;
      rows_per_block = 8;
    }
    else {
      threads_per_row = 16;
      rows_per_block = 16;
    }
    const size_type num_blocks =
      (row_count + rows_per_block -1 ) / rows_per_block;

    // Setup thread blocks and grid
    const dim3 dBlock( threads_per_row , rows_per_block , 1 );
    const dim3 dGrid( num_blocks, 1, 1 );

    // Setup shared memory for storing matrix values and column indices
    // Number of columns we process at a time -- making this bigger
    // requires more shared memory and reduces occupancy
    const int BlockSize = 32;
    if (sizeof(matrix_scalar) > 4)
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    else
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    const size_t shared =
      (BlockSize*dBlock.y) * (sizeof(size_type) + sizeof(matrix_scalar));

    // Launch kernel
    MeanFullOccupancyKernelLaunch<<<dGrid, dBlock, shared >>>
      ( Kernel<BlockSize>( A, x, y, a, b ) );
  }
};

/*
 * Disable this specialization as it is actually slower than the default
 * implementation of launching the single column kernel, one column at at time.
 * It looks like it uses a few more registers which is reducing occupancy.

// Kernel implementing y = A * x where PCE size of A is 1
//   A == Kokkos::CrsMatrix< Sacado::UQ::PCE<...>,...>, with A.values.sacado_size() == 1
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 2
template <typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class MeanMultiply< Kokkos::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                       MatrixOrdinal,
                                       Kokkos::Cuda,
                                       MatrixMemory,
                                       MatrixSize >,
                    Kokkos::View< Sacado::UQ::PCE<InputStorage>**,
                                  InputP... >,
                    Kokkos::View< Sacado::UQ::PCE<OutputStorage>**,
                                  OutputP... >,
                    >
{
public:
  typedef Sacado::UQ::PCE<MatrixStorage> MatrixValue;
  typedef Sacado::UQ::PCE<InputStorage> InputVectorValue;
  typedef Sacado::UQ::PCE<OutputStorage> OutputVectorValue;

  typedef Kokkos::Cuda MatrixDevice;
  typedef MatrixDevice execution_space;
  typedef Kokkos::CrsMatrix< const MatrixValue,
                             MatrixOrdinal,
                             MatrixDevice,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename MatrixValue::ordinal_type size_type;
  typedef Kokkos::View< InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;

  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename MatrixValue::value_type matrix_scalar;
  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  template <int BlockSize>
  struct Kernel {
    typedef Device execution_space;
    typedef typename Kokkos::FlatArrayType<matrix_values_type>::type matrix_array_type;
    typedef typename input_vector_type::array_type input_array_type;
    typedef typename output_vector_type::array_type output_array_type;

    const matrix_array_type   m_A_values ;
    const matrix_graph_type   m_A_graph ;
    const output_array_type   v_y ;
    const input_array_type    v_x ;
    const input_scalar        m_a ;
    const output_scalar       m_b ;
    const size_type           m_row_count;
    const size_type           m_num_col;
    const size_type           dim ;

    Kernel( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y ,
            const input_scalar & a ,
            const output_scalar & b )
      : m_A_values( A.values )
      , m_A_graph( A.graph )
      , v_y( y )
      , v_x( x )
      , m_a( a )
      , m_b( b )
      , m_row_count( A.graph.row_map.extent(0)-1 )
      , m_num_col( x.extent(1) )
      , dim( dimension_scalar(x) )
      {}

    __device__ void operator()(void) const
    {
      // Store matrix values and column indices in shared memory
      // to reduce global memory reads
      volatile matrix_scalar * const sh_A =
        kokkos_impl_cuda_shared_memory<matrix_scalar>();
      volatile size_type * const sh_col =
        reinterpret_cast<volatile size_type*>(sh_A + BlockSize*blockDim.y);

      // Check for valid row
      const size_type iBlockRow = blockDim.y*blockIdx.x + threadIdx.y;
      if (iBlockRow < m_row_count) {

        const size_type iEntryBegin = m_A_graph.row_map[ iBlockRow ];
        const size_type iEntryEnd   = m_A_graph.row_map[ iBlockRow + 1 ];

        // Initialize result
        if (m_b == output_scalar(0))
          for ( size_type xcol = 0; xcol < m_num_col; xcol++)
            for ( size_type pce = threadIdx.x; pce < dim ; pce+=blockDim.x )
              v_y(pce, iBlockRow, xcol) = 0.0;
        else
          for ( size_type xcol = 0; xcol < m_num_col; xcol++)
            for ( size_type pce = threadIdx.x; pce < dim ; pce+=blockDim.x )
              v_y(pce, iBlockRow, xcol) = m_b*v_y(pce, iBlockRow, xcol);

        // Loop over columns in chunks of size BlockSize
        for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
             col_block+=BlockSize) {
          const size_type num_col = col_block+BlockSize <= iEntryEnd ?
            BlockSize : iEntryEnd-col_block;

          // Read BlockSize entries column indices at a time to maintain
          // coalesced accesses
          for (size_type col=threadIdx.x; col<num_col; col+=blockDim.x) {
            sh_col[col*blockDim.y+threadIdx.y] =
              m_A_graph.entries(col_block+col);
            sh_A[col*blockDim.y+threadIdx.y] =
              m_A_values(col_block+col);
          }
          if (blockDim.x > Kokkos::Impl::CudaTraits::WarpSize)
            __syncthreads();

          // Do portion mat-vec for each PCE coefficient and for columns
          // within this block
          for ( size_type xcol = 0; xcol < m_num_col; xcol++) {
            for ( size_type pce = threadIdx.x; pce < dim ; pce+=blockDim.x ) {
              output_scalar s = 0.0;
              for ( size_type col = 0; col < num_col; ++col ) {
                const size_type iCol = sh_col[col*blockDim.y+threadIdx.y];
                const matrix_scalar aA = m_a*sh_A[col*blockDim.y+threadIdx.y];
                s += aA*v_x(pce, iCol, xcol);
              }
              v_y(pce, iBlockRow, xcol) += s;
            }
          }
        }
      }
    }
  };

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const size_t row_count = A.graph.row_map.extent(0) - 1;
    const size_type dim = dimension_scalar(x);

    // Compute number of threads for PCE coefficients and number of
    // matrix rows processed by each CUDA block.  A total of 256 threads
    // gives best occupancy
    size_type threads_per_row;
    size_type rows_per_block;
    if (dim >= 32) {
      threads_per_row = 32;
      rows_per_block = 6;
    }
    else {
      threads_per_row = 16;
      rows_per_block = 12;
    }
    const size_type num_blocks =
      (row_count + rows_per_block -1 ) / rows_per_block;

    // Setup thread blocks and grid
    const dim3 dBlock( threads_per_row , rows_per_block , 1 );
    const dim3 dGrid( num_blocks, 1, 1 );

    // Setup shared memory for storing matrix values and column indices
    // Number of columns we process at a time -- making this bigger
    // requires more shared memory and reduces occupancy
    const int BlockSize = 32;
    if (sizeof(matrix_scalar) > 4)
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    else
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    const size_t shared =
      (BlockSize*dBlock.y) * (sizeof(size_type) + sizeof(matrix_scalar));

    // Launch kernel
    // MeanFullOccupancyKernelLaunch<<<dGrid, dBlock, shared >>>
    //   ( Kernel<BlockSize>( A, x, y, a, b ) );

    Kokkos::Impl::cuda_parallel_launch_local_memory<<<dGrid, dBlock, shared >>>
      ( Kernel<BlockSize>( A, x, y, a, b ) );
  }
};
*/

} // namespace Stokhos

#endif /* #if defined( __CUDACC__) */

#endif /* #ifndef KOKKOS_CRSMATRIX_UQ_PCE_CUDA_HPP */
