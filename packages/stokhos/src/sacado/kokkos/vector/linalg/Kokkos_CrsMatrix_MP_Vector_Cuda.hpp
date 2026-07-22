// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_CRSMATRIX_MP_VECTOR_CUDA_HPP
#define KOKKOS_CRSMATRIX_MP_VECTOR_CUDA_HPP

#if defined( __CUDACC__)

#include "Kokkos_CrsMatrix_MP_Vector.hpp"
#include "Kokkos_Core.hpp"
#include "Stokhos_Cuda_DeviceProp.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::CrsMatrix for Sacado::MP::Vector scalar type
// and Cuda device
//----------------------------------------------------------------------------

namespace Stokhos {

namespace details {

// We can make things slightly faster for the ensemble multiply kernel with
// NumPerThread == 1 by using at most 32 registers per thread for 100%
// occupancy, which we get with these launch bounds on Kepler.  For
// NumPerThread > 1 it is worse though due to too many spilled registers.
template <typename Kernel>
__global__ void
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1024,2)
#endif
FullOccupancyKernelLaunch(Kernel kernel) {
  kernel();
}

// Kernel implementing y = A * x for Cuda device where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses the underlying 2-D view directly.
// Currently only works for statically sized MP::Vector
template <typename MatrixMemorySpace,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP,
          typename Update>
class MPMultiply< KokkosSparse::CrsMatrix<const Sacado::MP::Vector<MatrixStorage>,
                                    MatrixOrdinal,
                                    Kokkos::Device<Kokkos::Cuda, MatrixMemorySpace>,
                                    MatrixMemory,
                                    MatrixSize>,
                  Kokkos::View< const Sacado::MP::Vector<InputStorage>*,
                                InputP... >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                OutputP... >,
                  Update,
                  void
                >
{
public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Kokkos::Device<Kokkos::Cuda, MatrixMemorySpace> Device;
  typedef typename Device::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix<const MatrixValue,
                            MatrixOrdinal,
                            Device,
                            MatrixMemory,
                            MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;
  typedef Update update_type;

  // Multiply that uses shared memory for coalesced access of sparse
  // graph row offsets and column indices and Rank-1 vectors
  template <unsigned NumPerThread>
  struct Kernel {
    typedef typename output_vector_type::value_type output_vector_value;
    typedef typename output_vector_value::value_type scalar_type;

    const matrix_graph_type  m_Agraph;
    const matrix_values_type m_Avals;
    const input_vector_type  m_x;
    const output_vector_type m_y;
    const update_type m_update;
    const size_type m_row_count;

    Kernel( const matrix_type & A,
            const input_vector_type & x,
            const output_vector_type & y,
            const update_type& update )
      : m_Agraph( A.graph )
      , m_Avals( A.values )
      , m_x( x )
      , m_y( y )
      , m_update( update )
      , m_row_count( A.graph.row_map.extent(0)-1 )
      {}

    __device__
    inline void operator()(void) const
    {
      volatile size_type * const sh_row =
        kokkos_impl_cuda_shared_memory<size_type>();
      volatile size_type * const sh_col = sh_row + blockDim.y+1;

      const size_type tid = blockDim.x*threadIdx.y + threadIdx.x;
      const size_type nid = blockDim.x*blockDim.y;

      const size_type block_row = blockDim.y*blockIdx.x;

      // Read blockDim.y+1 row offsets in coalesced manner
      const size_type num_row =
        block_row+blockDim.y+1 <= m_row_count+1 ? blockDim.y+1 :
        m_row_count+1 - block_row;
      for (size_type i=tid; i<num_row; i+=nid)
        sh_row[i] = m_Agraph.row_map[block_row+i];
      __syncthreads();

      const size_type iRow = block_row + threadIdx.y;
      if (iRow < m_row_count) {
        scalar_type sum[NumPerThread];
        const size_type iEntryBegin = sh_row[threadIdx.y];
        const size_type iEntryEnd =   sh_row[threadIdx.y+1];

        for (size_type e=0; e<NumPerThread; ++e)
          sum[e] = 0;

        for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
             col_block+=blockDim.x) {
          const size_type num_col =
            col_block+blockDim.x <= iEntryEnd ?
              blockDim.x : iEntryEnd-col_block;

          // Read blockDim.x entries column indices at a time to maintain
          // coalesced accesses (don't need __syncthreads() assuming
          // blockDim.x <= warp_size
          // Note:  it might be a little faster if we ensured aligned access
          // to m_A.graph.entries() and m_A.values() below.
          if (threadIdx.x < num_col)
            sh_col[tid] = m_Agraph.entries(col_block+threadIdx.x);
          if (blockDim.x > Kokkos::Impl::CudaTraits::WarpSize)
            __syncthreads();

          for ( size_type col = 0; col < num_col; ++col ) {
            size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

            for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
                 ++e, ee+=blockDim.x) {
              sum[e] += m_Avals(col_block+col).fastAccessCoeff(ee) *
                        m_x(iCol).fastAccessCoeff(ee);
            }
          }

        }

        for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
             ++e, ee+=blockDim.x) {
          m_update( m_y(iRow).fastAccessCoeff(ee), sum[e] );
        }
      }
    } // operator()
  };

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    // Compute CUDA block and grid sizes.
    //
    // For Kepler, the best block size appears to be 256 threads with
    // 16 threads per vector for double precision, yielding 16 rows per
    // block.  Due to register usage, this gives 64 or 48 warps per SM
    // and thus 8 or 6 blocks per SM.  We use these values by default if
    // the user-specified block dimensions are zero

    const size_type value_dimension = Kokkos::dimension_scalar(x);

    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0)
      threads_per_vector = value_dimension ;
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0)
      rows_per_block = 256 / threads_per_vector;
    const size_type row_count = A.graph.row_map.extent(0)-1;
    const size_type num_blocks = (row_count+rows_per_block-1)/rows_per_block;
    const dim3 block( threads_per_vector, rows_per_block, 1 );
    const dim3 grid( num_blocks, 1 );

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = value_dimension / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != value_dimension, std::logic_error,
      "Entries/thread * threads/vector must equal number of vector entries");

    // Check threads_per_vector is not greater than warp size (kernels assume
    // this)
    const size_type warp_size = Kokkos::Impl::CudaTraits::WarpSize;
    TEUCHOS_TEST_FOR_EXCEPTION(
      threads_per_vector > warp_size, std::logic_error,
      "Threads/vector cannont exceed Cuda warp size");

    // Launch kernel based on static number of entries per thread
    if (num_per_thread == 1) {
      launch_impl<1>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 2) {
      launch_impl<2>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 3) {
      launch_impl<3>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 4) {
      launch_impl<4>( A, x, y, update, block, grid );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }

private:

  // Function to launch our kernel, templated on number of entries per thread
  // and shared memory choice
  template <unsigned num_per_thread>
  static void launch_impl( const matrix_type & A,
                           const input_vector_type & x,
                           const output_vector_type & y,
                           const update_type & update,
                           dim3 block,
                           dim3 grid)
  {
    typedef Kernel<num_per_thread> Krnl;

    // Use this to check occupancy, 64 is 100% on Kepler
    const bool occupancy_check = false;
    if (occupancy_check) {
      DeviceProp device_prop;
      size_type warps_per_sm;
      if (num_per_thread == 1)
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          FullOccupancyKernelLaunch<Krnl>);
      else
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          Kokkos::Impl::cuda_parallel_launch_local_memory<Krnl>);
      std::cout << "warps_per_sm = " << warps_per_sm
                << " max_warps_per_sm = " << device_prop.max_warps_per_sm
                << std::endl;
    }

    const size_t shared = (block.y+1 + block.x*block.y)*sizeof(size_type);
    if (sizeof(size_type) <= 4)
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    else
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    // Launch
    if (num_per_thread == 1)
      FullOccupancyKernelLaunch<<< grid, block, shared >>>
        ( Krnl( A, x, y, update ) );
    else
      Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
        ( Krnl( A, x, y, update ) );
  }
};

// Kernel implementing y = A * x for Cuda device where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses the underlying 2-D view directly.
// Currently only works for statically sized MP::Vector
template <typename MatrixMemorySpace,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP,
          typename Update>
class MPMultiply< KokkosSparse::CrsMatrix<const Sacado::MP::Vector<MatrixStorage>,
                                    MatrixOrdinal,
                                    Kokkos::Device<Kokkos::Cuda, MatrixMemorySpace>,
                                    MatrixMemory,
                                    MatrixSize>,
                  Kokkos::View< const Sacado::MP::Vector<InputStorage>**,
                                InputP... >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                OutputP... >,
                  Update,
                  void
                >
{
public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Kokkos::Device<Kokkos::Cuda, MatrixMemorySpace> Device;
  typedef typename Device::execution_space execution_space;
  typedef typename execution_space::size_type size_type;


  typedef KokkosSparse::CrsMatrix<const MatrixValue,
                            MatrixOrdinal,
                            Device,
                            MatrixMemory,
                            MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;
  typedef Update update_type;

  // Specialization that uses shared memory for coalesced access of sparse
  // graph row offsets and column indices and Rank-2 vectors
  template <unsigned NumPerThread>
  struct Kernel {
    typedef typename output_vector_type::value_type output_vector_value;
    typedef typename output_vector_value::value_type scalar_type;

    const matrix_graph_type  m_Agraph;
    const matrix_values_type m_Avals;
    const input_vector_type  m_x;
    const output_vector_type m_y;
    const update_type m_update;
    const size_type m_row_count;
    const size_type m_num_vec_cols;

    Kernel( const matrix_type & A,
            const input_vector_type & x,
            const output_vector_type & y,
            const update_type& update )
      : m_Agraph( A.graph )
      , m_Avals( A.values )
      , m_x( x )
      , m_y( y )
      , m_update( update )
      , m_row_count( A.graph.row_map.extent(0)-1 )
      , m_num_vec_cols( x.extent(1) )
      {}

    __device__
    inline void operator()(void) const
    {
      volatile size_type * const sh_row =
        kokkos_impl_cuda_shared_memory<size_type>();
      volatile size_type * const sh_col = sh_row + blockDim.y+1;

      const size_type tid = blockDim.x*threadIdx.y + threadIdx.x;
      const size_type nid = blockDim.x*blockDim.y;

      const size_type block_row = blockDim.y*blockIdx.x;

      // Read blockDim.y+1 row offsets in coalesced manner
      const size_type num_row =
        block_row+blockDim.y+1 <= m_row_count+1 ? blockDim.y+1 :
        m_row_count+1 - block_row;
      for (size_type i=tid; i<num_row; i+=nid)
        sh_row[i] = m_Agraph.row_map[block_row+i];
      __syncthreads();

      const size_type iRow = block_row + threadIdx.y;
      if (iRow < m_row_count) {
        scalar_type sum[NumPerThread];
        const size_type iEntryBegin = sh_row[threadIdx.y];
        const size_type iEntryEnd =   sh_row[threadIdx.y+1];

        // We could potentially improve performance by putting this loop
        // inside the matrix column loop, however that would require storing
        // an array for sum over vector columns, which would require shared
        // memory.
        for (size_type vec_col=0; vec_col<m_num_vec_cols; vec_col++) {

          for (size_type e=0; e<NumPerThread; ++e)
            sum[e] = 0;

          for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
               col_block+=blockDim.x) {
            const size_type num_col =
              col_block+blockDim.x <= iEntryEnd ?
                blockDim.x : iEntryEnd-col_block;

            // Read blockDim.x entries column indices at a time to maintain
            // coalesced accesses (don't need __syncthreads() assuming
            // blockDim.x <= warp_size
            // Note:  it might be a little faster if we ensured aligned access
            // to m_A.graph.entries() and m_A.values() below.
            if (threadIdx.x < num_col)
              sh_col[tid] = m_Agraph.entries(col_block+threadIdx.x);
            if (blockDim.x > Kokkos::Impl::CudaTraits::WarpSize)
              __syncthreads();

            for ( size_type col = 0; col < num_col; ++col ) {
              size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

              for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
                   ++e, ee+=blockDim.x) {
                sum[e] += m_Avals(col_block+col).fastAccessCoeff(ee) *
                          m_x(iCol, vec_col).fastAccessCoeff(ee);
              }
            }

          }

          for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
               ++e, ee+=blockDim.x) {
            m_update( m_y(iRow, vec_col).fastAccessCoeff(ee), sum[e] );
          }

        }
      }
    } // operator()
  };

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    // Compute CUDA block and grid sizes.
    //
    // For Kepler, the best block size appears to be 256 threads with
    // 16 threads per vector for double precision, yielding 16 rows per
    // block.  Due to register usage, this gives 64 or 48 warps per SM
    // and thus 8 or 6 blocks per SM.  We use these values by default if
    // the user-specified block dimensions are zero

    const size_type value_dimension = dimension_scalar(x);

    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0)
      threads_per_vector = value_dimension ;
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0)
      rows_per_block = 256 / threads_per_vector;
    const size_type row_count = A.graph.row_map.extent(0)-1;
    const size_type num_blocks = (row_count+rows_per_block-1)/rows_per_block;
    const dim3 block( threads_per_vector, rows_per_block, 1 );
    const dim3 grid( num_blocks, 1 );

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = value_dimension / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != value_dimension, std::logic_error,
      "Entries/thread * threads/vector must equal number of vector entries");

    // Launch kernel based on static number of entries per thread
    if (num_per_thread == 1) {
      launch_impl<1>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 2) {
      launch_impl<2>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 3) {
      launch_impl<3>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 4) {
      launch_impl<4>( A, x, y, update, block, grid );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }

private:

  // Function to launch our kernel, templated on number of entries per thread
  // and shared memory choice
  template <unsigned num_per_thread>
  static void launch_impl( const matrix_type & A,
                           const input_vector_type & x,
                           const output_vector_type & y,
                           const update_type & update,
                           dim3 block,
                           dim3 grid)
  {
    typedef Kernel<num_per_thread> Krnl;

    // Use this to check occupancy, 64 is 100% on Kepler
    const bool occupancy_check = false;
    if (occupancy_check) {
      DeviceProp device_prop;
      size_type warps_per_sm;
      if (num_per_thread == 1)
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          FullOccupancyKernelLaunch<Krnl>);
      else
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          Kokkos::Impl::cuda_parallel_launch_local_memory<Krnl>);
      std::cout << "warps_per_sm = " << warps_per_sm
                << " max_warps_per_sm = " << device_prop.max_warps_per_sm
                << std::endl;
    }

    const size_t shared = (block.y+1 + block.x*block.y)*sizeof(size_type);
    if (sizeof(size_type) <= 4)
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    else
      cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    // Launch
    if (num_per_thread == 1)
      FullOccupancyKernelLaunch<<< grid, block, shared >>>
        ( Krnl( A, x, y, update ) );
    else
      Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
        ( Krnl( A, x, y, update ) );
  }
};

} // namespace details

} // namespace Stokhos

#endif /* #if defined( __CUDACC__) */

#endif /* #ifndef KOKKOS_CRSMATRIX_MP_VECTOR_CUDA_HPP */
