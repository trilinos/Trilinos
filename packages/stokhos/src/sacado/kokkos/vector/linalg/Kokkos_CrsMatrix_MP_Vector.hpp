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
#include "Kokkos_CrsMatrix.hpp"

#include "Kokkos_Parallel.hpp"
#include "Stokhos_Multiply.hpp"

// For computing ParallelWorkRequest
#include "Kokkos_hwloc.hpp"
#include "Kokkos_Cuda.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::CrsMatrix for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

//! Specialization of SparseRowView<> for Sacado::MP::Vector scalar type
/*!
 * Here we store the view and offset directly instead of pointers, since
 * pointers don't currently work well with View< Sacado::MP::Vector<...>,... >
 */
template <typename Storage,
          typename OrdinalType,
          typename MemoryTraits,
          typename SizeType>
struct SparseRowView< Kokkos::CrsMatrix< Sacado::MP::Vector<Storage>,
                                         OrdinalType,
                                         typename Storage::device_type,
                                         MemoryTraits,
                                         SizeType> > {
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<Storage>,
                             OrdinalType,
                             typename Storage::device_type,
                             MemoryTraits,
                             SizeType> MatrixType;
  typedef typename MatrixType::values_type values_type;
  typedef typename MatrixType::index_type index_type;
  typedef typename values_type::sacado_mp_vector_view_type scalar_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

private:

  values_type values_;
  index_type colidx_;
  const int offset_, stride_;

public:

  KOKKOS_INLINE_FUNCTION
  SparseRowView (const values_type& values,
                 const index_type& colidx,
                 const int stride,
                 const int count,
                 const int offset) :
    values_ (values), colidx_ (colidx),
    offset_ (offset), stride_ (stride), length (count)
  {}

  const int length;

  KOKKOS_INLINE_FUNCTION
  scalar_type & value (const int& i) const {
    return values_(offset_+i*stride_);
  }

  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx (const int& i) const {
    return colidx_[offset_+i*stride_];
  }
};

//! Specialization of SparseRowViewConst<> for Sacado::MP::Vector scalar type
/*!
 * Here we store the view and offset directly instead of pointers, since
 * pointers don't currently work well with View< Sacado::MP::Vector<...>,... >
 */
template <typename Storage,
          typename OrdinalType,
          typename MemoryTraits,
          typename SizeType>
struct SparseRowViewConst< Kokkos::CrsMatrix< Sacado::MP::Vector<Storage>,
                                              OrdinalType,
                                              typename Storage::device_type,
                                              MemoryTraits,
                                              SizeType> > {
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<Storage>,
                             OrdinalType,
                             typename Storage::device_type,
                             MemoryTraits,
                             SizeType> MatrixType;
  typedef typename MatrixType::values_type values_type;
  typedef typename MatrixType::index_type index_type;
  typedef typename values_type::sacado_mp_vector_view_type scalar_type;
  typedef const typename MatrixType::non_const_ordinal_type ordinal_type;

private:

  values_type values_;
  index_type colidx_;
  const int offset_, stride_;

public:

  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst (const values_type& values,
                      const index_type& colidx,
                      const int stride,
                      const int count,
                      const int offset) :
    values_ (values), colidx_ (colidx),
    offset_ (offset), stride_ (stride), length (count)
  {}

  const int length;

  KOKKOS_INLINE_FUNCTION
  scalar_type& value (const int& i) const {
    return values_(offset_+i*stride_);
  }

  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx (const int& i) const {
    return colidx_[offset_+i*stride_];
  }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
// Sparse matrix-vector product kernels
//----------------------------------------------------------------------------

namespace Stokhos {

// Specialization of ViewRank<> for View< Sacado::MP::Vector<...>,...> since
// the rank enum for the view specialization isn't right
template <typename Storage, typename L, typename D, typename M, typename S>
struct ViewRank< Kokkos::View<Sacado::MP::Vector<Storage>*,L,D,M,S> > {
  typedef IntegralRank<1> type;
};

namespace details {

/*
 * Compute work range = (begin, end) such that adjacent threads/blocks write to
 * separate cache lines
 */
template <typename scalar_type, typename device_type, typename size_type>
KOKKOS_INLINE_FUNCTION
size_type
compute_work_range( const size_type work_count,
                    const size_type thread_count,
                    const size_type thread_rank,
                    size_type& work_begin)
{
  enum { cache_line =
         Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value ? 128 : 64 };
  enum { work_align = cache_line / sizeof(scalar_type) };
  enum { work_shift = Kokkos::Impl::power_of_two< work_align >::value };
  enum { work_mask  = work_align - 1 };

  const size_type work_per_thread =
    ( ( ( ( work_count + work_mask ) >> work_shift ) + thread_count - 1 ) /
      thread_count ) << work_shift ;

  work_begin = thread_rank * work_per_thread;
  size_type work_end = work_begin + work_per_thread;
  if (work_begin > work_count)
    work_begin = work_count;
  if (work_end > work_count)
    work_end = work_count;
  return work_end;
}

} // namespace details

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputMemory,
          typename OutputStorage,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                              Kokkos::LayoutRight,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                              Kokkos::LayoutRight,
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
                        Kokkos::LayoutRight,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        Kokkos::LayoutRight,
                        Device,
                        OutputMemory > output_vector_type;

private:

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
    output_vector_type  m_y;

    MultiplyKernel( const matrix_type & A,
                    const input_vector_type & x,
                    output_vector_type & y )
      : m_A( A )
      , m_x( x )
      , m_y( y )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( device_type dev ) const
    {
      // 2-D distribution of threads: num_vector_threads x num_row_threads
      // where the x-dimension are vector threads and the y dimension are
      // row threads
      const size_type num_vector_threads = m_A.dev_config.block_dim.x ;
      const size_type num_row_threads    = m_A.dev_config.block_dim.y ;
      const size_type row_rank = dev.team_rank() / num_vector_threads ;

      // Create local views with corresponding offset into the vector
      // dimension based on vector_rank and reduced number of vector entries

      // Partition Sacado::MP::Vector as [ NumPerThread * rank .. NumPerThread * ( rank + 1 ) )
      const Sacado::MP::VectorPartition part( NumPerThread * dev.team_rank() ,
                                              NumPerThread * ( dev.team_rank() + 1 ) );

      const matrix_local_view_type A = Kokkos::subview<matrix_local_view_type>( m_A.values, part );
      const input_local_view_type  x = Kokkos::subview<input_local_view_type>(  m_x , part );
      const output_local_view_type y = Kokkos::subview<output_local_view_type>( m_y , part );

      // const matrix_values_type A(m_A.values);
      // const input_vector_type x(m_x);
      // const output_vector_type y(m_y);

      // Compute range of rows processed for each thread block
      const size_type row_count = m_A.graph.row_map.dimension_0()-1;
      size_type work_begin;
      const size_type work_end =
        details::compute_work_range<typename scalar_type::value_type,device_type,size_type>(row_count, dev.league_size(), dev.league_rank(), work_begin);

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
      for (size_type iBlockRow=work_begin+row_rank*rows_per_thread;
           iBlockRow<work_end; iBlockRow+=row_block_size) {

        // Loop over rows within block
        const size_type row_end = iBlockRow+rows_per_thread <= work_end ?
          rows_per_thread : work_end - iBlockRow;
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
          y(iRow) = sum;

        } // row loop

      } // block row loop

    } // operator()

  }; // MatrixKernel<NumPerThread>

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     output_vector_type & y )
  {
    const bool is_cuda =
      Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value;

    // By default, use one one entry per thread for GPU and
    // one thread per vector for CPU/MIC
    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0) {
      if (is_cuda)
        threads_per_vector = x.dimension_1();
      else
        threads_per_vector = 1;
    }

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = x.dimension_1() / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != x.dimension_1(), std::logic_error,
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
      Kokkos::parallel_for( config, MultiplyKernel<1>(A,x,y) );
    else if (num_per_thread == 2)
      Kokkos::parallel_for( config, MultiplyKernel<2>(A,x,y) );
    else if (num_per_thread == 3)
      Kokkos::parallel_for( config, MultiplyKernel<3>(A,x,y) );
    else if (num_per_thread == 4)
      Kokkos::parallel_for( config, MultiplyKernel<4>(A,x,y) );
    else if (num_per_thread == 8)
      Kokkos::parallel_for( config, MultiplyKernel<8>(A,x,y) );
    else if (num_per_thread == 12)
      Kokkos::parallel_for( config, MultiplyKernel<12>(A,x,y) );
    else if (num_per_thread == 16)
      Kokkos::parallel_for( config, MultiplyKernel<16>(A,x,y) );
    else if (num_per_thread == 32)
      Kokkos::parallel_for( config, MultiplyKernel<32>(A,x,y) );
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
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
          typename InputMemory,
          typename OutputStorage,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                              Kokkos::LayoutRight,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                              Kokkos::LayoutRight,
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
                        Kokkos::LayoutRight,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        Kokkos::LayoutRight,
                        Device,
                        OutputMemory > output_vector_type;

private:

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
    typename output_vector_type::array_type  m_y;

    MultiplyKernel( const matrix_type & A,
                    const input_vector_type & x,
                    output_vector_type & y )
      : m_A( A )
      , m_Avals( A.values )
      , m_x( x )
      , m_y( y )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( device_type dev ) const
    {
      // 2-D distribution of threads: num_vector_threads x num_row_threads
      // where the x-dimension are vector threads and the y dimension are
      // row threads
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
      size_type work_begin;
      const size_type work_end =
        details::compute_work_range<scalar_type,device_type,size_type>(
          row_count, dev.league_size(), dev.league_rank(), work_begin);
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
      for (size_type iBlockRow=work_begin+row_rank*rows_per_thread;
           iBlockRow<work_end; iBlockRow+=row_block_size) {

        // Loop over rows within block
        const size_type row_end = iBlockRow+rows_per_thread <= work_end ?
          rows_per_thread : work_end - iBlockRow;
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
            y[e*stride] = sum[e];

        } // row loop

      } // block row loop

    } // operator()

  }; // MatrixKernel<NumPerThread>

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     output_vector_type & y )
  {
    const bool is_cuda =
      Kokkos::Impl::is_same<device_type,Kokkos::Cuda>::value;

    // By default, use one one entry per thread for GPU and
    // one thread per vector for CPU/MIC
    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0) {
      if (is_cuda)
        threads_per_vector = x.dimension_1();
      else
        threads_per_vector = 1;
    }

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = x.dimension_1() / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != x.dimension_1(), std::logic_error,
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
      Kokkos::parallel_for( config, MultiplyKernel<1>(A,x,y) );
    else if (num_per_thread == 2)
      Kokkos::parallel_for( config, MultiplyKernel<2>(A,x,y) );
    else if (num_per_thread == 3)
      Kokkos::parallel_for( config, MultiplyKernel<3>(A,x,y) );
    else if (num_per_thread == 4)
      Kokkos::parallel_for( config, MultiplyKernel<4>(A,x,y) );
    else if (num_per_thread == 8)
      Kokkos::parallel_for( config, MultiplyKernel<8>(A,x,y) );
    else if (num_per_thread == 12)
      Kokkos::parallel_for( config, MultiplyKernel<12>(A,x,y) );
    else if (num_per_thread == 16)
      Kokkos::parallel_for( config, MultiplyKernel<16>(A,x,y) );
    else if (num_per_thread == 32)
      Kokkos::parallel_for( config, MultiplyKernel<32>(A,x,y) );
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
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
  Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
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
  //typedef Stokhos::DefaultMultiply Tag;
  typedef Stokhos::EnsembleMultiply Tag;
  Stokhos::multiply(A, x, y, Tag());
}

}

#endif /* #ifndef KOKKOS_CRSMATRIX_MP_VECTOR_HPP */
