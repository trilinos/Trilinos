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

#ifndef KOKKOS_CRSMATRIX_MP_VECTOR_HPP
#define KOKKOS_CRSMATRIX_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_MV_MP_Vector.hpp" // for some utilities

#include "Kokkos_Parallel.hpp"
#include "Stokhos_Multiply.hpp"

// For computing ParallelWorkRequest
#include "Kokkos_hwloc.hpp"
#include "Kokkos_Cuda.hpp"

#include "Teuchos_TestForException.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::CrsMatrix for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Stokhos {

namespace details {

template <typename Matrix, typename InputVector, typename OutputVector,
          typename Update = MultiplyAssign>
class MPMultiply {};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename Update>
class MPMultiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                     MatrixOrdinal,
                                     Device,
                                     MatrixMemory,
                                     MatrixSize>,
                  Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                                InputLayout,
                                Device,
                                InputMemory >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                OutputLayout,
                                Device,
                                OutputMemory >,
                  Update
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue*,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;
  typedef Update update_type;

  template <unsigned NumPerThread>
  struct MultiplyKernel {

    typedef Device device_type;

    // Typenames for thread-local views
    typedef typename Kokkos::LocalMPVectorView<matrix_values_type,
                                               NumPerThread>::type matrix_local_view_type;
    typedef typename Kokkos::LocalMPVectorView<input_vector_type,
                                               NumPerThread>::type input_local_view_type;
    typedef typename Kokkos::LocalMPVectorView<output_vector_type,
                                               NumPerThread>::type output_local_view_type;

    typedef typename output_local_view_type::value_type scalar_type;

    const matrix_type  m_A;
    const input_vector_type  m_x;
    const output_vector_type  m_y;
    const update_type m_update;

    MultiplyKernel( const matrix_type & A,
                    const input_vector_type & x,
                    const output_vector_type & y,
                    const update_type& update )
      : m_A( A )
      , m_x( x )
      , m_y( y )
      , m_update( update )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( device_type dev ) const
    {
      // 2-D distribution of threads: num_vector_threads x num_row_threads
      // where the x-dimension are vector threads and the y dimension are
      // row threads
      //
      // Note:  The actual team size granted by Kokkos may not match the
      // request, which we need to handle in the kernel launch
      // (because this kernel is effectively templated on the team size).
      // Currently Kokkos doesn't make this easy.
      const size_type num_vector_threads = m_A.dev_config.block_dim.x ;
      const size_type num_row_threads    = m_A.dev_config.block_dim.y ;
      const size_type row_rank = dev.team_rank() / num_vector_threads ;
      const size_type vec_rank = dev.team_rank() % num_vector_threads ;

#if !defined(__CUDA_ARCH__)
      TEUCHOS_TEST_FOR_EXCEPTION(
        num_vector_threads * num_row_threads != size_type(dev.team_size()),
        std::logic_error,
        "num_vector_threads (" << num_vector_threads <<
        ") * num_row_threads (" << num_row_threads <<
        ") != dev.team_size (" << dev.team_size() << ")!");
      TEUCHOS_TEST_FOR_EXCEPTION(
        m_A.dev_config.num_blocks != size_type(dev.league_size()),
        std::logic_error,
        "num_blocks (" << m_A.dev_config.num_blocks <<
        ") != dev.league_size (" << dev.league_size() << ")!");
#endif

      // Create local views with corresponding offset into the vector
      // dimension based on vector_rank and reduced number of vector entries

      // Partition Sacado::MP::Vector as
      // [ NumPerThread * rank .. NumPerThread * ( rank + 1 ) )
      const Sacado::MP::VectorPartition part( NumPerThread * vec_rank ,
                                              NumPerThread * (vec_rank + 1 ) );

      const matrix_local_view_type A =
        Kokkos::subview<matrix_local_view_type>( m_A.values, part );
      const input_local_view_type  x =
        Kokkos::subview<input_local_view_type>(  m_x , part );
      const output_local_view_type y =
        Kokkos::subview<output_local_view_type>( m_y , part );

      // const matrix_values_type A(m_A.values);
      // const input_vector_type x(m_x);
      // const output_vector_type y(m_y);

      // Compute range of rows processed for each thread block
      const size_type row_count = m_A.graph.row_map.dimension_0()-1;
      const size_type league_size = dev.league_size();
      const size_type league_rank = dev.league_rank();
      const Kokkos::pair<size_type,size_type> work_range =
        details::compute_work_range<typename scalar_type::value_type>(
          dev, row_count, league_size, league_rank);

      // To make better use of L1 cache on the CPU/MIC, we move through the
      // row range where each thread processes a cache-line's worth of rows,
      // with adjacent threads processing adjacent cache-lines.
      // For Cuda, adjacent threads process adjacent rows.
      const size_type cache_line =
        Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value ? 1 : 64;
      const size_type scalar_size = sizeof(scalar_type);
      const size_type rows_per_thread = (cache_line+scalar_size-1)/scalar_size;
      const size_type row_block_size = rows_per_thread * num_row_threads;

      scalar_type sum;

      // Loop over rows in blocks of row_block_size
      for (size_type iBlockRow=work_range.first+row_rank*rows_per_thread;
           iBlockRow<work_range.second; iBlockRow+=row_block_size) {

        // Loop over rows within block
        const size_type row_end =
          iBlockRow+rows_per_thread <= work_range.second ?
            rows_per_thread : work_range.second - iBlockRow;
        for (size_type row=0; row<row_end; ++row) {
          const size_type iRow = iBlockRow + row;

          // Compute mat-vec for this row
          const size_type iEntryBegin = m_A.graph.row_map[iRow];
          const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];
          sum = 0.0;
          for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
            size_type iCol = m_A.graph.entries(iEntry);
            sum += A(iEntry) * x(iCol);
          }
          m_update( y(iRow), sum );

        } // row loop

      } // block row loop

    } // operator()

  }; // MatrixKernel<NumPerThread>

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const bool is_cuda =
      Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value;

    matrix_type AA = A;

    // By default, use one one entry per thread for GPU and
    // one thread per vector for CPU/MIC
    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0) {
      if (is_cuda)
        threads_per_vector = x.sacado_size();
      else
        threads_per_vector = 1;
    }
    AA.dev_config.block_dim.x = threads_per_vector;

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = x.sacado_size() / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != x.sacado_size(),
      std::logic_error,
      "Entries/thread * threads/vector must equal number of vector entries");

    // By default, use a block size of 256 for GPU and number of hyperthreads
    // per core for CPU/MIC
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0) {
      if (is_cuda)
        rows_per_block = 256 / threads_per_vector;
      else {
        rows_per_block =
          Kokkos::hwloc::get_available_threads_per_core();
      }
    }
    AA.dev_config.block_dim.y = rows_per_block;

    // Thread team size
    size_type team_size = threads_per_vector * rows_per_block;

    // Number of teams -- For GPU, each block does rows_per_block rows,
    // for CPU/MIC, league_size is the number of cores
    size_type league_size = A.dev_config.num_blocks;
    if (league_size == 0) {
      if (is_cuda) {
        const size_type row_count = A.graph.row_map.dimension_0()-1;
        league_size = (row_count+rows_per_block-1)/rows_per_block;
      }
      else
        league_size =
          Kokkos::hwloc::get_available_numa_count() *
          Kokkos::hwloc::get_available_cores_per_numa();
    }
    AA.dev_config.num_blocks = league_size;

    // Parallel launch with corresponding number of vector entries per thread
    Kokkos::ParallelWorkRequest config(league_size, team_size);
    if (num_per_thread == 1)
      Kokkos::parallel_for( config, MultiplyKernel<1>(AA,x,y,update) );
    else if (num_per_thread == 2)
      Kokkos::parallel_for( config, MultiplyKernel<2>(AA,x,y,update) );
    else if (num_per_thread == 3)
      Kokkos::parallel_for( config, MultiplyKernel<3>(AA,x,y,update) );
    else if (num_per_thread == 4)
      Kokkos::parallel_for( config, MultiplyKernel<4>(AA,x,y,update) );
    else if (num_per_thread == 8)
      Kokkos::parallel_for( config, MultiplyKernel<8>(AA,x,y,update) );
    else if (num_per_thread == 12)
      Kokkos::parallel_for( config, MultiplyKernel<12>(AA,x,y,update) );
    else if (num_per_thread == 16)
      Kokkos::parallel_for( config, MultiplyKernel<16>(AA,x,y,update) );
    else if (num_per_thread == 20)
      Kokkos::parallel_for( config, MultiplyKernel<20>(AA,x,y,update) );
    else if (num_per_thread == 24)
      Kokkos::parallel_for( config, MultiplyKernel<24>(AA,x,y,update) );
    else if (num_per_thread == 32)
      Kokkos::parallel_for( config, MultiplyKernel<32>(AA,x,y,update) );
    else if (num_per_thread == 40)
      Kokkos::parallel_for( config, MultiplyKernel<40>(AA,x,y,update) );
    else if (num_per_thread == 48)
      Kokkos::parallel_for( config, MultiplyKernel<48>(AA,x,y,update) );
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }
};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename Update>
class MPMultiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                     MatrixOrdinal,
                                     Device,
                                     MatrixMemory,
                                     MatrixSize>,
                  Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                                InputLayout,
                                Device,
                                InputMemory >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                OutputLayout,
                                Device,
                                OutputMemory >,
                  Update
                  >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue**,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;
  typedef Update update_type;

  template <unsigned NumPerThread>
  struct MultiplyKernel {

    typedef Device device_type;

    // Typenames for thread-local views
    typedef typename Kokkos::LocalMPVectorView<matrix_values_type,
                                               NumPerThread>::type matrix_local_view_type;
    typedef typename Kokkos::LocalMPVectorView<input_vector_type,
                                               NumPerThread>::type input_local_view_type;
    typedef typename Kokkos::LocalMPVectorView<output_vector_type,
                                               NumPerThread>::type output_local_view_type;

    typedef typename output_local_view_type::value_type scalar_type;

    const matrix_type  m_A;
    const input_vector_type  m_x;
    const output_vector_type  m_y;
    const update_type m_update;

    MultiplyKernel( const matrix_type & A,
                    const input_vector_type & x,
                    const output_vector_type & y,
                    const update_type& update )
      : m_A( A )
      , m_x( x )
      , m_y( y )
      , m_update( update )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( device_type dev ) const
    {
      // 2-D distribution of threads: num_vector_threads x num_row_threads
      // where the x-dimension are vector threads and the y dimension are
      // row threads
      //
      // Note:  The actual team size granted by Kokkos may not match the
      // request, which we need to handle in the kernel launch
      // (because this kernel is effectively templated on the team size).
      // Currently Kokkos doesn't make this easy.
      const size_type num_vector_threads = m_A.dev_config.block_dim.x ;
      const size_type num_row_threads    = m_A.dev_config.block_dim.y ;
      const size_type row_rank = dev.team_rank() / num_vector_threads ;
      const size_type vec_rank = dev.team_rank() % num_vector_threads ;

#if !defined(__CUDA_ARCH__)
      TEUCHOS_TEST_FOR_EXCEPTION(
        num_vector_threads * num_row_threads != size_type(dev.team_size()),
        std::logic_error,
        "num_vector_threads (" << num_vector_threads <<
        ") * num_row_threads (" << num_row_threads <<
        ") != dev.team_size (" << dev.team_size() << ")!");
      TEUCHOS_TEST_FOR_EXCEPTION(
        m_A.dev_config.num_blocks != size_type(dev.league_size()),
        std::logic_error,
        "num_blocks (" << m_A.dev_config.num_blocks <<
        ") != dev.league_size (" << dev.league_size() << ")!");
#endif

      // Create local views with corresponding offset into the vector
      // dimension based on vector_rank and reduced number of vector entries

      // Partition Sacado::MP::Vector as
      // [ NumPerThread * rank .. NumPerThread * ( rank + 1 ) )
      const Sacado::MP::VectorPartition part( NumPerThread * vec_rank ,
                                              NumPerThread * (vec_rank + 1 ) );

      const matrix_local_view_type A =
        Kokkos::subview<matrix_local_view_type>( m_A.values, part );
      const input_local_view_type  x =
        Kokkos::subview<input_local_view_type>(  m_x , part );
      const output_local_view_type y =
        Kokkos::subview<output_local_view_type>( m_y , part );

      // const matrix_values_type A(m_A.values);
      // const input_vector_type x(m_x);
      // const output_vector_type y(m_y);

      // Compute range of rows processed for each thread block
      const size_type row_count = m_A.graph.row_map.dimension_0()-1;
      const size_type league_size = dev.league_size();
      const size_type league_rank = dev.league_rank();
      const Kokkos::pair<size_type,size_type> work_range =
        details::compute_work_range<typename scalar_type::value_type>(
          dev, row_count, league_size, league_rank);

      // To make better use of L1 cache on the CPU/MIC, we move through the
      // row range where each thread processes a cache-line's worth of rows,
      // with adjacent threads processing adjacent cache-lines.
      // For Cuda, adjacent threads process adjacent rows.
      const size_type cache_line =
        Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value ? 1 : 64;
      const size_type scalar_size = sizeof(scalar_type);
      const size_type rows_per_thread = (cache_line+scalar_size-1)/scalar_size;
      const size_type row_block_size = rows_per_thread * num_row_threads;

      const size_type num_col = m_y.dimension_1();

      scalar_type sum;

      // Loop over rows in blocks of row_block_size
      for (size_type iBlockRow=work_range.first+row_rank*rows_per_thread;
           iBlockRow< work_range.second; iBlockRow+=row_block_size) {

        // Loop over rows within block
        const size_type row_end =
          iBlockRow+rows_per_thread <=  work_range.second ?
            rows_per_thread :  work_range.second - iBlockRow;
        for (size_type row=0; row<row_end; ++row) {
          const size_type iRow = iBlockRow + row;

          // Loop over columns of x, y
          for (size_type col=0; col<num_col; ++col) {

            // Compute mat-vec for this row
            const size_type iEntryBegin = m_A.graph.row_map[iRow];
            const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];
            sum = 0.0;
            for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
              size_type iCol = m_A.graph.entries(iEntry);
              sum += A(iEntry) * x(iCol,col);
            }
            m_update( y(iRow,col), sum );

          } // x, y column loop

        } // row loop

      } // block row loop

    } // operator()

  }; // MatrixKernel<NumPerThread>

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const bool is_cuda =
      Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value;

    matrix_type AA = A;

    // By default, use one one entry per thread for GPU and
    // one thread per vector for CPU/MIC
    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0) {
      if (is_cuda)
        threads_per_vector = x.sacado_size();
      else
        threads_per_vector = 1;
    }
    AA.dev_config.block_dim.x = threads_per_vector;

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = x.sacado_size() / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != x.sacado_size(),
      std::logic_error,
      "Entries/thread * threads/vector must equal number of vector entries");

    // By default, use a block size of 256 for GPU and number of hyperthreads
    // per core for CPU/MIC
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0) {
      if (is_cuda)
        rows_per_block = 256 / threads_per_vector;
      else {
        rows_per_block =
          Kokkos::hwloc::get_available_threads_per_core();
      }
    }
    AA.dev_config.block_dim.y = rows_per_block;

    // Thread team size
    size_type team_size = threads_per_vector * rows_per_block;

    // Number of teams -- For GPU, each block does rows_per_block rows,
    // for CPU/MIC, league_size is the number of cores
    size_type league_size = A.dev_config.num_blocks;
    if (league_size == 0) {
      if (is_cuda) {
        const size_type row_count = A.graph.row_map.dimension_0()-1;
        league_size = (row_count+rows_per_block-1)/rows_per_block;
      }
      else
        league_size =
          Kokkos::hwloc::get_available_numa_count() *
          Kokkos::hwloc::get_available_cores_per_numa();
    }
    AA.dev_config.num_blocks = league_size;

    // Parallel launch with corresponding number of vector entries per thread
    Kokkos::ParallelWorkRequest config(league_size, team_size);
    if (num_per_thread == 1)
      Kokkos::parallel_for( config, MultiplyKernel<1>(AA,x,y,update) );
    else if (num_per_thread == 2)
      Kokkos::parallel_for( config, MultiplyKernel<2>(AA,x,y,update) );
    else if (num_per_thread == 3)
      Kokkos::parallel_for( config, MultiplyKernel<3>(AA,x,y,update) );
    else if (num_per_thread == 4)
      Kokkos::parallel_for( config, MultiplyKernel<4>(AA,x,y,update) );
    else if (num_per_thread == 8)
      Kokkos::parallel_for( config, MultiplyKernel<8>(AA,x,y,update) );
    else if (num_per_thread == 12)
      Kokkos::parallel_for( config, MultiplyKernel<12>(AA,x,y,update) );
    else if (num_per_thread == 16)
      Kokkos::parallel_for( config, MultiplyKernel<16>(AA,x,y,update) );
    else if (num_per_thread == 20)
      Kokkos::parallel_for( config, MultiplyKernel<20>(AA,x,y,update) );
    else if (num_per_thread == 24)
      Kokkos::parallel_for( config, MultiplyKernel<24>(AA,x,y,update) );
    else if (num_per_thread == 32)
      Kokkos::parallel_for( config, MultiplyKernel<32>(AA,x,y,update) );
    else if (num_per_thread == 40)
      Kokkos::parallel_for( config, MultiplyKernel<40>(AA,x,y,update) );
    else if (num_per_thread == 48)
      Kokkos::parallel_for( config, MultiplyKernel<48>(AA,x,y,update) );
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }
};

template <typename Matrix, typename InputVector, typename OutputVector,
          typename Update = MultiplyAssign>
class EnsembleMultiply {};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses the underlying 2-D view directly.
// Currently only works for statically sized MP::Vector
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename Update>
class EnsembleMultiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                           MatrixOrdinal,
                                           Device,
                                           MatrixMemory,
                                           MatrixSize>,
                        Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                                      InputLayout,
                                      Device,
                                      InputMemory >,
                        Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                      OutputLayout,
                                      Device,
                                      OutputMemory >,
                        Update
                      >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue*,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;
  typedef Update update_type;

  template <unsigned NumPerThread>
  struct MultiplyKernel {

    typedef Device device_type;

    typedef typename OutputVectorValue::value_type scalar_type;
    typedef typename InputVectorValue::value_type x_scalar_type;
    typedef typename OutputVectorValue::value_type y_scalar_type;
    typedef typename MatrixValue::value_type A_scalar_type;

    const matrix_type  m_A;
    const typename matrix_values_type::array_type m_Avals;
    const typename input_vector_type::array_type  m_x;
    const typename output_vector_type::array_type  m_y;
    const update_type m_update;

    MultiplyKernel( const matrix_type & A,
                    const input_vector_type & x,
                    const output_vector_type & y,
                    const update_type& update )
      : m_A( A )
      , m_Avals( A.values )
      , m_x( x )
      , m_y( y )
      , m_update( update )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( device_type dev ) const
    {
      // 2-D distribution of threads: num_vector_threads x num_row_threads
      // where the x-dimension are vector threads and the y dimension are
      // row threads
      //
      // Note:  The actual team size granted by Kokkos may not match the
      // request, which we need to handle in the kernel launch
      // (because this kernel is effectively templated on the team size).
      // Currently Kokkos doesn't make this easy.
      const size_type num_vector_threads = m_A.dev_config.block_dim.x;
      const size_type num_row_threads = m_A.dev_config.block_dim.y;
      const size_type vector_rank = dev.team_rank() % num_vector_threads;
      const size_type row_rank = dev.team_rank() / num_vector_threads;

      // We have to extract pointers to the A, x, and y views below
      // because the Intel compiler does not appear to be able to vectorize
      // through them.  Thus we need to compute the correct stride for Cuda.
      using Kokkos::Impl::is_same;
      const bool is_cuda = is_same<device_type, Kokkos::Cuda>::value;
      const size_type stride = is_cuda ? num_vector_threads : size_type(1);

      // Compute range of rows processed for each thread block
      const size_type row_count = m_A.graph.row_map.dimension_0()-1;
      const size_type league_size = dev.league_size();
      const size_type league_rank = dev.league_rank();
      const Kokkos::pair<size_type,size_type> work_range =
        details::compute_work_range<scalar_type>(
          dev, row_count, league_size, league_rank);
      const size_type vector_offset =
        is_cuda ? vector_rank : vector_rank*NumPerThread;

      // To make better use of L1 cache on the CPU/MIC, we move through the
      // row range where each thread processes a cache-line's worth of rows,
      // with adjacent threads processing adjacent cache-lines.
      // For Cuda, adjacent threads process adjacent rows.
      const size_type cache_line = is_cuda ? 1 : 64;
      const size_type scalar_size = sizeof(scalar_type);
      const size_type rows_per_thread = (cache_line+scalar_size-1)/scalar_size;
      const size_type row_block_size = rows_per_thread * num_row_threads;

      scalar_type sum[NumPerThread];

      // Loop over rows in blocks of row_block_size
      for (size_type iBlockRow=work_range.first+row_rank*rows_per_thread;
           iBlockRow<work_range.second; iBlockRow+=row_block_size) {

        // Loop over rows within block
        const size_type row_end =
          iBlockRow+rows_per_thread <= work_range.second ?
            rows_per_thread : work_range.second - iBlockRow;
        for (size_type row=0; row<row_end; ++row) {
          const size_type iRow = iBlockRow + row;

          // Compute mat-vec for this row
          const size_type iEntryBegin = m_A.graph.row_map[iRow];
          const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

          y_scalar_type * const y = &m_y(iRow,vector_offset);

          for (size_type e=0; e<NumPerThread; ++e)
            sum[e] = 0;

          for ( size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry ) {
            size_type iCol = m_A.graph.entries(iEntry);

            const A_scalar_type * const A = &m_Avals(iEntry,vector_offset);
            const x_scalar_type * const x = &m_x(iCol,vector_offset);

            for (size_type e=0; e<NumPerThread; ++e)
              sum[e] += A[e*stride] * x[e*stride];
          }

          for (size_type e=0; e<NumPerThread; ++e)
            m_update( y[e*stride], sum[e] );

        } // row loop

      } // block row loop

    } // operator()

  }; // MatrixKernel<NumPerThread>

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const bool is_cuda =
      Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value;

    // By default, use one one entry per thread for GPU and
    // one thread per vector for CPU/MIC
    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0) {
      if (is_cuda)
        threads_per_vector = x.sacado_size();
      else
        threads_per_vector = 1;
    }

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = x.sacado_size() / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != x.sacado_size(),
      std::logic_error,
      "Entries/thread * threads/vector must equal number of vector entries");

    // By default, use a block size of 256 for GPU and number of hyperthreads
    // per core for CPU/MIC
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0) {
      if (is_cuda)
        rows_per_block = 256 / threads_per_vector;
      else {
        rows_per_block =
          Kokkos::hwloc::get_available_threads_per_core();
      }
    }

    // Thread team size
    size_type team_size = threads_per_vector * rows_per_block;

    // Number of teams -- For GPU, each block does rows_per_block rows,
    // for CPU/MIC, league_size is the number of cores
    size_type league_size = A.dev_config.num_blocks;
    if (league_size == 0) {
      if (is_cuda) {
        const size_type row_count = A.graph.row_map.dimension_0()-1;
        league_size = (row_count+rows_per_block-1)/rows_per_block;
      }
      else
        league_size =
          Kokkos::hwloc::get_available_numa_count() *
          Kokkos::hwloc::get_available_cores_per_numa();
    }

    // Parallel launch with corresponding number of vector entries per thread
    Kokkos::ParallelWorkRequest config(league_size, team_size);
    if (num_per_thread == 1)
      Kokkos::parallel_for( config, MultiplyKernel<1>(A,x,y,update) );
    else if (num_per_thread == 2)
      Kokkos::parallel_for( config, MultiplyKernel<2>(A,x,y,update) );
    else if (num_per_thread == 3)
      Kokkos::parallel_for( config, MultiplyKernel<3>(A,x,y,update) );
    else if (num_per_thread == 4)
      Kokkos::parallel_for( config, MultiplyKernel<4>(A,x,y,update) );
    else if (num_per_thread == 8)
      Kokkos::parallel_for( config, MultiplyKernel<8>(A,x,y,update) );
    else if (num_per_thread == 12)
      Kokkos::parallel_for( config, MultiplyKernel<12>(A,x,y,update) );
    else if (num_per_thread == 16)
      Kokkos::parallel_for( config, MultiplyKernel<16>(A,x,y,update) );
    else if (num_per_thread == 20)
      Kokkos::parallel_for( config, MultiplyKernel<20>(A,x,y,update) );
    else if (num_per_thread == 24)
      Kokkos::parallel_for( config, MultiplyKernel<24>(A,x,y,update) );
    else if (num_per_thread == 32)
      Kokkos::parallel_for( config, MultiplyKernel<32>(A,x,y,update) );
    else if (num_per_thread == 40)
      Kokkos::parallel_for( config, MultiplyKernel<40>(A,x,y,update) );
    else if (num_per_thread == 48)
      Kokkos::parallel_for( config, MultiplyKernel<48>(A,x,y,update) );
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }
};

template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename UpdateType = MultiplyAssign>
struct MultiplyType {
  typedef MPMultiply<MatrixType,InputVectorType,OutputVectorType,UpdateType> type;
};

} // namespace details

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                              InputLayout,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                              OutputLayout,
                              Device,
                              OutputMemory >
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue*,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y )
  {
    typedef details::MPMultiply<matrix_type,input_vector_type,output_vector_type> multiply_type;
    multiply_type::apply(A,x,y,details::MultiplyAssign());
  }
};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                              InputLayout,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                              OutputLayout,
                              Device,
                              OutputMemory >
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue**,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y )
  {
    typedef details::MPMultiply<matrix_type,input_vector_type,output_vector_type> multiply_type;
    multiply_type::apply(A,x,y,details::MultiplyAssign());
  }
};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses the underlying 2-D view directly.
// Currently only works for statically sized MP::Vector
struct EnsembleMultiply {};
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                              InputLayout,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                              OutputLayout,
                              Device,
                              OutputMemory >,
                void,
                IntegralRank<1>,
                EnsembleMultiply >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue*,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y )
  {
    typedef details::EnsembleMultiply<matrix_type,input_vector_type,output_vector_type> multiply_type;
    multiply_type::apply(A,x,y,details::MultiplyAssign());
  }
};

// Overload of multiply() for ensemble tag
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType>
void multiply(const MatrixType& A,
              const InputVectorType& x,
              OutputVectorType& y,
              EnsembleMultiply tag) {
  typedef Multiply<MatrixType,InputVectorType,OutputVectorType,void,typename ViewRank<InputVectorType>::type,EnsembleMultiply> multiply_type;
  multiply_type::apply( A, x, y );
}

} // namespace Stokhos

namespace Kokkos {

// Overload of Kokkos::MV_Multiply for Sacado::MP::Vector scalar types
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
                                          OutputLayout,
                                          Device,
                                          OutputMemory >& y,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;
  typedef typename Stokhos::details::MultiplyType<MatrixType,InputVectorType,
    OutputVectorType>::type multiply_type;

  multiply_type::apply( A, x, y, Stokhos::details::MultiplyAssign() );
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
                     OutputLayout,
                     Device,
                     OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef typename InputStorage::value_type value_type;
  typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;

  if (!Impl::is_mp_vector_constant(a)) {
    Impl::raise_mp_vector_error(
      "MV_Multiply not implemented for non-constant a");
  }

  value_type aa = a.fastAccessCoeff(0);
  if (aa == value_type(1)) {
    // y = A*x
    typedef Stokhos::details::MultiplyAssign UpdateType;
    typedef typename Stokhos::details::MultiplyType<MatrixType,InputVectorType,
      OutputVectorType,UpdateType>::type multiply_type;
    multiply_type::apply( A, x, y, UpdateType() );
  }
  else {
    // y = a*A*x
    typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
    typedef typename Stokhos::details::MultiplyType<MatrixType,InputVectorType,
      OutputVectorType,UpdateType>::type multiply_type;
    multiply_type::apply( A, x, y, UpdateType(aa) );
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Sacado::MP::Vector<InputStorage>& b,
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
                     OutputLayout,
                     Device,
                     OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef typename InputStorage::value_type value_type;
  typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;

  if (!Impl::is_mp_vector_constant(a) || !Impl::is_mp_vector_constant(b)) {
    Impl::raise_mp_vector_error(
      "MV_Multiply not implemented for non-constant a or b");
  }

  value_type aa = a.fastAccessCoeff(0);
  value_type bb = b.fastAccessCoeff(0);
  if (bb == value_type(0)) {
    if (aa == value_type(1)) {
      // y = A*x
      typedef Stokhos::details::MultiplyAssign UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y = a*A*x
      typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
  else if (bb == value_type(1)) {
    if (aa == value_type(1)) {
      // y += A*x
      typedef Stokhos::details::MultiplyUpdate UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y += a*A*x
      typedef Stokhos::details::MultiplyScaledUpdate<value_type> UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
  else {
    // y = a*A*x + b*y
    typedef Stokhos::details::MultiplyScaledUpdate2<value_type> UpdateType;
    typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
    multiply_type::apply( A, x, y, UpdateType(aa,bb) );
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(y_1D, A, x_1D);
  }
  else {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
      OutputLayout, Device, OutputMemory > OutputVectorType;
    typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
      MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
      InputLayout, Device, InputMemory > InputVectorType;

    // y = A*x
    typedef Stokhos::details::MultiplyAssign UpdateType;
    typedef typename Stokhos::details::MultiplyType<MatrixType,
      InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
    multiply_type::apply( A, x, y, UpdateType() );
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(y_1D, a, A, x_1D);
  }
  else {
    typedef typename InputStorage::value_type value_type;
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
      OutputLayout, Device, OutputMemory > OutputVectorType;
    typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
      MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
      InputLayout, Device, InputMemory > InputVectorType;

    if (!Impl::is_mp_vector_constant(a)) {
      Impl::raise_mp_vector_error(
        "MV_Multiply not implemented for non-constant a");
    }

    value_type aa = a.fastAccessCoeff(0);
    if (aa == value_type(1)) {
      // y = A*x
      typedef Stokhos::details::MultiplyAssign UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y = a*A*x
      typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Sacado::MP::Vector<InputStorage>& b,
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(b, y_1D, a, A, x_1D);
  }
  else {
    typedef typename InputStorage::value_type value_type;
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
      OutputLayout, Device, OutputMemory > OutputVectorType;
    typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
      MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
      InputLayout, Device, InputMemory > InputVectorType;

    if (!Impl::is_mp_vector_constant(a) || !Impl::is_mp_vector_constant(b)) {
      Impl::raise_mp_vector_error(
        "MV_Multiply not implemented for non-constant a or b");
    }

    value_type aa = a.fastAccessCoeff(0);
    value_type bb = b.fastAccessCoeff(0);
    if (bb == value_type(0)) {
      if (aa == value_type(1)) {
        // y = A*x
        typedef Stokhos::details::MultiplyAssign UpdateType;
        typedef typename Stokhos::details::MultiplyType<MatrixType,
          InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
        multiply_type::apply( A, x, y, UpdateType() );
      }
      else {
        // y = a*A*x
        typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
        typedef typename Stokhos::details::MultiplyType<MatrixType,
          InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
        multiply_type::apply( A, x, y, UpdateType(aa) );
      }
    }
    else if (bb == value_type(1)) {
      if (aa == value_type(1)) {
        // y += A*x
        typedef Stokhos::details::MultiplyUpdate UpdateType;
        typedef typename Stokhos::details::MultiplyType<MatrixType,
          InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
        multiply_type::apply( A, x, y, UpdateType() );
      }
      else {
        // y += a*A*x
        typedef Stokhos::details::MultiplyScaledUpdate<value_type> UpdateType;
        typedef typename Stokhos::details::MultiplyType<MatrixType,
          InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
        multiply_type::apply( A, x, y, UpdateType(aa) );
      }
    }
    else {
      // y = a*A*x + b*y
      typedef Stokhos::details::MultiplyScaledUpdate2<value_type> UpdateType;
      typedef typename Stokhos::details::MultiplyType<MatrixType,
        InputVectorType,OutputVectorType,UpdateType>::type multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa,bb) );
    }
  }
}

}

#endif /* #ifndef KOKKOS_CRSMATRIX_MP_VECTOR_HPP */
