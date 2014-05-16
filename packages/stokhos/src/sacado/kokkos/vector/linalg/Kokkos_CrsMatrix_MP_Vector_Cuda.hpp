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

#ifndef KOKKOS_CRSMATRIX_MP_VECTOR_CUDA_HPP
#define KOKKOS_CRSMATRIX_MP_VECTOR_CUDA_HPP

#if defined( __CUDACC__)

#include "Kokkos_CrsMatrix_MP_Vector.hpp"
#include "Kokkos_Cuda.hpp"
#include "Cuda/Kokkos_Cuda_Parallel.hpp"
#include "Stokhos_Cuda_DeviceProp.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::CrsMatrix for Sacado::MP::Vector scalar type
// and Cuda device
//----------------------------------------------------------------------------

namespace Stokhos {

namespace details {

// Kernels for ensemble matrix-vector multiply
// You would like these to be nested within the Multiply<> specializations
// below, but you can't partial specialize within class-scope.
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          unsigned NumPerThread,
          bool UseShared,
          unsigned Rank = InputVectorType::Rank,
          typename Layout = typename InputVectorType::array_layout,
          typename Update = MultiplyAssign>
struct EnsembleMultiplyKernel {};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices
// Rank-1 and LayoutRight vectors (note matrix is always LayoutRight)
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          unsigned NumPerThread,
          typename Update>
struct EnsembleMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              NumPerThread, true, 1, Kokkos::LayoutRight,
                              Update> {
  typedef MatrixType matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;
  typedef typename OutputVectorType::value_type OutputVectorValue;
  typedef typename OutputVectorValue::value_type scalar_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  const matrix_graph_type  m_Agraph;
  const typename matrix_values_type::array_type m_Avals;
  const typename input_vector_type::array_type  m_x;
  const typename output_vector_type::array_type m_y;
  const update_type m_update;
  const size_type m_row_count;

  EnsembleMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type& update )
    : m_Agraph( A.graph )
    , m_Avals( A.values )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
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
          col_block+blockDim.x <= iEntryEnd ? blockDim.x : iEntryEnd-col_block;

        // Read blockDim.x entries column indices at a time to maintain
        // coalesced accesses (don't need __syncthreads() assuming
        // blockDim.x <= warp_size
        // Note:  it might be a little faster if we ensured aligned access
        // to m_A.graph.entries() and m_A.values() below.
        if (threadIdx.x < num_col)
          sh_col[tid] = m_Agraph.entries(col_block+threadIdx.x);

        for ( size_type col = 0; col < num_col; ++col ) {
          size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

          for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
               ++e, ee+=blockDim.x) {
            sum[e] += m_Avals(col_block+col, ee) * m_x(iCol, ee);
          }
        }

      }

      for (size_type e=0, ee=threadIdx.x; e<NumPerThread; ++e, ee+=blockDim.x) {
        m_update( m_y(iRow, ee), sum[e] );
      }
    }
  } // operator()
};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices
// Rank-1 and LayoutLeft vectors (note matrix is always LayoutRight)
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          unsigned NumPerThread,
          typename Update>
struct EnsembleMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              NumPerThread, true, 1, Kokkos::LayoutLeft,
                              Update> {
  typedef MatrixType matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;
  typedef typename OutputVectorType::value_type OutputVectorValue;
  typedef typename OutputVectorValue::value_type scalar_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  const matrix_graph_type  m_Agraph;
  const typename matrix_values_type::array_type m_Avals;
  const typename input_vector_type::array_type  m_x;
  const typename output_vector_type::array_type m_y;
  const update_type m_update;
  const size_type m_row_count;

  EnsembleMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type& update )
    : m_Agraph( A.graph )
    , m_Avals( A.values )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
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
          col_block+blockDim.x <= iEntryEnd ? blockDim.x : iEntryEnd-col_block;

        // Read blockDim.x entries column indices at a time to maintain
        // coalesced accesses (don't need __syncthreads() assuming
        // blockDim.x <= warp_size
        // Note:  it might be a little faster if we ensured aligned access
        // to m_A.graph.entries() and m_A.values() below.
        if (threadIdx.x < num_col)
          sh_col[tid] = m_Agraph.entries(col_block+threadIdx.x);

        for ( size_type col = 0; col < num_col; ++col ) {
          size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

          for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
               ++e, ee+=blockDim.x) {
            sum[e] += m_Avals(col_block+col, ee) * m_x(ee, iCol);
          }
        }

      }

      for (size_type e=0, ee=threadIdx.x; e<NumPerThread; ++e, ee+=blockDim.x) {
        m_update( m_y(ee, iRow), sum[e] );
      }
    }
  } // operator()
};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices
// Rank-2 and LayoutLeft vectors (note matrix is always LayoutRight)
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          unsigned NumPerThread,
          typename Update>
struct EnsembleMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              NumPerThread, true, 2, Kokkos::LayoutLeft,
                              Update> {
  typedef MatrixType matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;
  typedef typename OutputVectorType::value_type OutputVectorValue;
  typedef typename OutputVectorValue::value_type scalar_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  const matrix_graph_type  m_Agraph;
  const typename matrix_values_type::array_type m_Avals;
  const typename input_vector_type::array_type  m_x;
  const typename output_vector_type::array_type m_y;
  const update_type m_update;
  const size_type m_row_count;
  const size_type m_num_vec_cols;

  EnsembleMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type& update )
    : m_Agraph( A.graph )
    , m_Avals( A.values )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
    , m_num_vec_cols( x.dimension_1() )
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

          for ( size_type col = 0; col < num_col; ++col ) {
            size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

            for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
                 ++e, ee+=blockDim.x) {
              sum[e] += m_Avals(col_block+col, ee) * m_x(ee, iCol, vec_col);
            }
          }

        }

        for (size_type e=0, ee=threadIdx.x; e<NumPerThread;
             ++e, ee+=blockDim.x) {
          m_update( m_y(ee, iRow, vec_col), sum[e] );
        }

      }
    }
  } // operator()
};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices and a single ensemble value per thread
// Rank-1 and LayoutRight vectors (note matrix is always LayoutRight)
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename Update>
struct EnsembleMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              1, true, 1, Kokkos::LayoutRight,
                              Update> {
  typedef MatrixType matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;
  typedef typename OutputVectorType::value_type OutputVectorValue;
  typedef typename OutputVectorValue::value_type scalar_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  const matrix_graph_type  m_Agraph;
  const typename matrix_values_type::array_type m_Avals;
  const typename input_vector_type::array_type  m_x;
  const typename output_vector_type::array_type m_y;
  const update_type m_update;
  const size_type m_row_count;

  EnsembleMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type & update )
    : m_Agraph( A.graph )
    , m_Avals( A.values )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
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
      const size_type iEntryBegin = sh_row[threadIdx.y];
      const size_type iEntryEnd =   sh_row[threadIdx.y+1];

      scalar_type sum = 0;

      for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
           col_block+=blockDim.x) {
        const size_type num_col =
          col_block+blockDim.x <= iEntryEnd ? blockDim.x : iEntryEnd-col_block;

        // Read blockDim.x entries column indices at a time to maintain
        // coalesced accesses (don't need __syncthreads() assuming
        // blockDim.x <= warp_size
        // Note:  it might be a little faster if we ensured aligned access
        // to m_A.graph.entries() and m_A.values() below.
        if (threadIdx.x < num_col)
          sh_col[tid] = m_Agraph.entries(col_block+threadIdx.x);

        for ( size_type col = 0; col < num_col; ++col ) {
          size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

          sum += m_Avals(col_block+col, threadIdx.x) * m_x(iCol, threadIdx.x);
        }

      }

      m_update( m_y(iRow, threadIdx.x), sum );
    }
  } // operator()
};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices and a single ensemble value per thread
// Rank-1 and LayoutLeft vectors (note matrix is always LayoutRight)
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename Update>
struct EnsembleMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              1, true, 1, Kokkos::LayoutLeft,
                              Update> {
  typedef MatrixType matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;
  typedef typename OutputVectorType::value_type OutputVectorValue;
  typedef typename OutputVectorValue::value_type scalar_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  const matrix_graph_type  m_Agraph;
  const typename matrix_values_type::array_type m_Avals;
  const typename input_vector_type::array_type  m_x;
  const typename output_vector_type::array_type m_y;
  const update_type m_update;
  const size_type m_row_count;

  EnsembleMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type& update )
    : m_Agraph( A.graph )
    , m_Avals( A.values )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
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
      const size_type iEntryBegin = sh_row[threadIdx.y];
      const size_type iEntryEnd =   sh_row[threadIdx.y+1];

      scalar_type sum = 0;

      for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
           col_block+=blockDim.x) {
        const size_type num_col =
          col_block+blockDim.x <= iEntryEnd ? blockDim.x : iEntryEnd-col_block;

        // Read blockDim.x entries column indices at a time to maintain
        // coalesced accesses (don't need __syncthreads() assuming
        // blockDim.x <= warp_size
        // Note:  it might be a little faster if we ensured aligned access
        // to m_A.graph.entries() and m_A.values() below.
        if (threadIdx.x < num_col)
          sh_col[tid] = m_Agraph.entries(col_block+threadIdx.x);

        for ( size_type col = 0; col < num_col; ++col ) {
          size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

          sum += m_Avals(col_block+col, threadIdx.x) * m_x(threadIdx.x, iCol);
        }

      }

      m_update( m_y(threadIdx.x, iRow), sum );
    }
  } // operator()
};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices and a single ensemble value per thread
// Rank-2 and LayoutLeft vectors (note matrix is always LayoutRight)
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename Update>
struct EnsembleMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              1, true, 2, Kokkos::LayoutLeft,
                              Update> {
  typedef MatrixType matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;
  typedef typename OutputVectorType::value_type OutputVectorValue;
  typedef typename OutputVectorValue::value_type scalar_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  const matrix_graph_type  m_Agraph;
  const typename matrix_values_type::array_type m_Avals;
  const typename input_vector_type::array_type  m_x;
  const typename output_vector_type::array_type m_y;
  const update_type m_update;
  const size_type m_row_count;
  const size_type m_num_vec_cols;

  EnsembleMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type& update )
    : m_Agraph( A.graph )
    , m_Avals( A.values )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
    , m_num_vec_cols( x.dimension_1() )
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
      const size_type iEntryBegin = sh_row[threadIdx.y];
      const size_type iEntryEnd =   sh_row[threadIdx.y+1];

      // We could potentially improve performance by putting this loop
      // inside the matrix column loop, however that would require storing
      // an array for sum over vector columns, which would require shared
      // memory.
      for (size_type vec_col=0; vec_col<m_num_vec_cols; vec_col++) {

        scalar_type sum = 0;

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

          for ( size_type col = 0; col < num_col; ++col ) {
            size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

            sum += m_Avals(col_block+col, threadIdx.x) *
              m_x(threadIdx.x, iCol, vec_col);
          }

        }

        m_update( m_y(threadIdx.x, iRow, vec_col), sum );

      }
    }
  } // operator()
};

/////////////////////////////////////////////////////

// Kernels for Sacado::MP::Vector matrix-vector multiply
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          unsigned NumPerThread,
          bool UseShared,
          unsigned Rank = InputVectorType::Rank,
          typename Update = MultiplyAssign>
struct MPVectorMultiplyKernel {};

// Specialization that uses shared memory for coalesced access of sparse
// graph row offsets and column indices
template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          unsigned NumPerThread,
          typename Update>
struct MPVectorMultiplyKernel<MatrixType,
                              InputVectorType,
                              OutputVectorType,
                              NumPerThread,
                              true, 1, Update> {
  typedef MatrixType matrix_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;

  typedef typename Kokkos::Cuda device_type;
  typedef typename device_type::size_type size_type;
  typedef Update update_type;

  typedef typename matrix_type::values_type matrix_values_type;

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
  const output_vector_type m_y;
  const update_type m_update;
  const size_type m_row_count;

  MPVectorMultiplyKernel( const matrix_type & A,
                          const input_vector_type & x,
                          const output_vector_type & y,
                          const update_type& update )
    : m_A( A )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    , m_row_count( A.graph.row_map.dimension_0()-1 )
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
      sh_row[i] = m_A.graph.row_map[block_row+i];
    __syncthreads();

    const size_type iRow = block_row + threadIdx.y;
    if (iRow < m_row_count) {

      const Sacado::MP::VectorPartition part( NumPerThread * threadIdx.x ,
                                              NumPerThread * ( threadIdx.x + 1 ) );

      const matrix_local_view_type A =
        Kokkos::subview<matrix_local_view_type>( m_A.values, part );
      const input_local_view_type  x =
        Kokkos::subview<input_local_view_type>(  m_x, part );
      const output_local_view_type y =
        Kokkos::subview<output_local_view_type>( m_y, part );

      const size_type iEntryBegin = sh_row[threadIdx.y];
      const size_type iEntryEnd =   sh_row[threadIdx.y+1];

      scalar_type sum = 0.0;

      for (size_type col_block=iEntryBegin; col_block<iEntryEnd;
           col_block+=blockDim.x) {
        const size_type num_col =
          col_block+blockDim.x <= iEntryEnd ? blockDim.x : iEntryEnd-col_block;

        // Read blockDim.x entries column indices at a time to maintain
        // coalesced accesses (don't need __syncthreads() assuming
        // blockDim.x <= warp_size
        // Note:  it might be a little faster if we ensured aligned access
        // to m_A.graph.entries() and m_A.values() below.
        if (threadIdx.x < num_col)
          sh_col[tid] = m_A.graph.entries(col_block+threadIdx.x);

        for ( size_type col = 0; col < num_col; ++col ) {
          size_type iCol = sh_col[blockDim.x*threadIdx.y + col];

          sum += A(col_block+col) * x(iCol);
        }

      }

      m_update( y(iRow), sum );
    }
  } // operator()
};

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
template <typename MatrixStorage,
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
class EnsembleMultiply< Kokkos::CrsMatrix<Sacado::MP::Vector<MatrixStorage>,
                                          MatrixOrdinal,
                                          Kokkos::Cuda,
                                          MatrixMemory,
                                          MatrixSize>,
                        Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                                      InputLayout,
                                      Kokkos::Cuda,
                                      InputMemory >,
                        Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                      OutputLayout,
                                      Kokkos::Cuda,
                                      OutputMemory >,
                        Update
                      >
{
public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef typename Kokkos::Cuda Device;
  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix<MatrixValue,
                            MatrixOrdinal,
                            device_type,
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

    const size_type value_dimension = x.sacado_size();

    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0)
      threads_per_vector = value_dimension ;
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0)
      rows_per_block = 256 / threads_per_vector;
    const size_type row_count = A.graph.row_map.dimension_0()-1;
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

    // The shared memory kernel is slightly faster
    const bool use_shared = true;

    // Launch kernel based on static number of entries per thread
    if (num_per_thread == 1) {
      launch_impl<1,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 2) {
      launch_impl<2,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 3) {
      launch_impl<3,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 4) {
      launch_impl<4,use_shared>( A, x, y, update, block, grid );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }

private:

  // Function to launch our kernel, templated on number of entries per thread
  // and shared memory choice
  template <unsigned num_per_thread, bool use_shared>
  static void launch_impl( const matrix_type & A,
                           const input_vector_type & x,
                           const output_vector_type & y,
                           const update_type & update,
                           dim3 block,
                           dim3 grid)
  {
    // The kernel we will launch
    typedef EnsembleMultiplyKernel<matrix_type,input_vector_type,output_vector_type,num_per_thread,use_shared,1,InputLayout,update_type> Kernel;

    // Use this to check occupancy, 64 is 100% on Kepler
    const bool occupancy_check = false;
    if (occupancy_check) {
      DeviceProp device_prop;
      size_type warps_per_sm;
      if (num_per_thread == 1)
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          FullOccupancyKernelLaunch<Kernel>);
      else
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          Kokkos::Impl::cuda_parallel_launch_local_memory<Kernel>);
      std::cout << "warps_per_sm = " << warps_per_sm
                << " max_warps_per_sm = " << device_prop.max_warps_per_sm
                << std::endl;
    }

    size_t shared = 0;
    if (use_shared) {
      shared = (block.y+1 + block.x*block.y)*sizeof(size_type);
      if (sizeof(size_type) <= 4)
        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
      else
        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    }

    // Launch
    if (num_per_thread == 1)
      FullOccupancyKernelLaunch<<< grid, block, shared >>>
        ( Kernel( A, x, y, update ) );
    else
      Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
        ( Kernel( A, x, y, update ) );
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
template <typename MatrixStorage,
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
class EnsembleMultiply< Kokkos::CrsMatrix<Sacado::MP::Vector<MatrixStorage>,
                                          MatrixOrdinal,
                                          Kokkos::Cuda,
                                          MatrixMemory,
                                          MatrixSize>,
                        Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                                      InputLayout,
                                      Kokkos::Cuda,
                                      InputMemory >,
                        Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                      OutputLayout,
                                      Kokkos::Cuda,
                                      OutputMemory >,
                        Update
                      >
{
public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef typename Kokkos::Cuda Device;
  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix<MatrixValue,
                            MatrixOrdinal,
                            device_type,
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

    const size_type value_dimension = x.sacado_size();

    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0)
      threads_per_vector = value_dimension ;
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0)
      rows_per_block = 256 / threads_per_vector;
    const size_type row_count = A.graph.row_map.dimension_0()-1;
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

    // The shared memory kernel is slightly faster
    const bool use_shared = true;

    // Launch kernel based on static number of entries per thread
    if (num_per_thread == 1) {
      launch_impl<1,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 2) {
      launch_impl<2,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 3) {
      launch_impl<3,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 4) {
      launch_impl<4,use_shared>( A, x, y, update, block, grid );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }

private:

  // Function to launch our kernel, templated on number of entries per thread
  // and shared memory choice
  template <unsigned num_per_thread, bool use_shared>
  static void launch_impl( const matrix_type & A,
                           const input_vector_type & x,
                           const output_vector_type & y,
                           const update_type & update,
                           dim3 block,
                           dim3 grid)
  {
    // The kernel we will launch
    typedef EnsembleMultiplyKernel<matrix_type,input_vector_type,output_vector_type,num_per_thread,use_shared,2,InputLayout,update_type> Kernel;

    // Use this to check occupancy, 64 is 100% on Kepler
    const bool occupancy_check = false;
    if (occupancy_check) {
      DeviceProp device_prop;
      size_type warps_per_sm;
      if (num_per_thread == 1)
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          FullOccupancyKernelLaunch<Kernel>);
      else
        warps_per_sm = device_prop.get_resident_warps_per_sm(
          Kokkos::Impl::cuda_parallel_launch_local_memory<Kernel>);
      std::cout << "warps_per_sm = " << warps_per_sm
                << " max_warps_per_sm = " << device_prop.max_warps_per_sm
                << std::endl;
    }

    size_t shared = 0;
    if (use_shared) {
      shared = (block.y+1 + block.x*block.y)*sizeof(size_type);
      if (sizeof(size_type) <= 4)
        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
      else
        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    }

    // Launch
    if (num_per_thread == 1)
      FullOccupancyKernelLaunch<<< grid, block, shared >>>
        ( Kernel( A, x, y, update ) );
    else
      Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
        ( Kernel( A, x, y, update ) );
  }
};

// Kernel implementing y = A * x for Cuda device where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixStorage,
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
class MPMultiply< Kokkos::CrsMatrix<Sacado::MP::Vector<MatrixStorage>,
                                    MatrixOrdinal,
                                    Kokkos::Cuda,
                                    MatrixMemory,
                                    MatrixSize>,
                  Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                                InputLayout,
                                Kokkos::Cuda,
                                InputMemory >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                OutputLayout,
                                Kokkos::Cuda,
                                OutputMemory >,
                  Update
                >
{
public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef typename Kokkos::Cuda Device;
  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix<MatrixValue,
                            MatrixOrdinal,
                            device_type,
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

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type& update )
  {
    // Compute CUDA block and grid sizes.
    //
    // For Kepler, the best block size appears to be 256 threads with
    // 16 threads per vector for double precision, yielding 16 rows per
    // block.  Due to register usage, this gives 64 or 48 warps per SM
    // and thus 8 or 6 blocks per SM.  We use these values by default if
    // the user-specified block dimensions are zero

    const size_type value_dimension = x.sacado_size();

    size_type threads_per_vector = A.dev_config.block_dim.x;
    if (threads_per_vector == 0)
      threads_per_vector = value_dimension ;
    size_type rows_per_block = A.dev_config.block_dim.y;
    if (rows_per_block == 0)
      rows_per_block = 256 / threads_per_vector;
    const size_type row_count = A.graph.row_map.dimension_0()-1;
    const size_type num_blocks = (row_count+rows_per_block-1)/rows_per_block;
    const dim3 block( threads_per_vector, rows_per_block, 1 );
    const dim3 grid( num_blocks, 1 );

    // Check threads_per_vector evenly divides number of vector entries
    size_type num_per_thread = value_dimension / threads_per_vector;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_per_thread * threads_per_vector != value_dimension, std::logic_error,
      "Entries/thread * threads/vector must equal number of vector entries");

    // The shared memory kernel is slightly faster
    const bool use_shared = true;

    // Launch kernel based on static number of entries per thread
    if (num_per_thread == 1) {
      launch_impl<1,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 2) {
      launch_impl<2,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 3) {
      launch_impl<3,use_shared>( A, x, y, update, block, grid );
    }
    else if (num_per_thread == 4) {
      launch_impl<4,use_shared>( A, x, y, update, block, grid );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid num_per_thread == " << num_per_thread);
  }

private:

  // Function to launch our kernel, templated on number of entries per thread
  // and shared memory choice
  template <unsigned num_per_thread, bool use_shared>
  static void launch_impl( const matrix_type & A,
                           const input_vector_type & x,
                           const output_vector_type & y,
                           const update_type& update,
                           dim3 block,
                           dim3 grid)
  {
    // The kernel we will launch
    typedef MPVectorMultiplyKernel<matrix_type,input_vector_type,output_vector_type,num_per_thread,use_shared,1,update_type> Kernel;

    // Use this to check occupancy, 64 is 100% on Kepler
    const bool occupancy_check = false;
    if (occupancy_check) {
      DeviceProp device_prop;
      size_type warps_per_sm = device_prop.get_resident_warps_per_sm(
          Kokkos::Impl::cuda_parallel_launch_local_memory<Kernel>);
      std::cout << "warps_per_sm = " << warps_per_sm
                << " max_warps_per_sm = " << device_prop.max_warps_per_sm
                << std::endl;
    }

    size_t shared = 0;
    if (use_shared) {
      shared = (block.y+1 + block.x*block.y)*sizeof(size_type);
      if (sizeof(size_type) <= 4)
        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
      else
        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    }

    // Launch
    Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
      ( Kernel( A, x, y, update ) );
  }
};

template <typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename UpdateType>
struct MultiplyType< Kokkos::CrsMatrix<Sacado::MP::Vector<MatrixStorage>,
                                          MatrixOrdinal,
                                          Kokkos::Cuda,
                                          MatrixMemory,
                                          MatrixSize>,
                        Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                                      InputLayout,
                                      Kokkos::Cuda,
                                      InputMemory >,
                        Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                      OutputLayout,
                                      Kokkos::Cuda,
                                      OutputMemory >,
                        UpdateType
                      > {
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                             MatrixOrdinal,
                             Kokkos::Cuda,
                             MatrixMemory,
                             MatrixSize > MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                        InputLayout,
                        Kokkos::Cuda,
                        InputMemory > InputVectorType;
  typedef Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                        OutputLayout,
                        Kokkos::Cuda,
                        OutputMemory > OutputVectorType;
  typedef EnsembleMultiply<MatrixType,InputVectorType,OutputVectorType,UpdateType> type;
};

template <typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename UpdateType>
struct MultiplyType< Kokkos::CrsMatrix<Sacado::MP::Vector<MatrixStorage>,
                                          MatrixOrdinal,
                                          Kokkos::Cuda,
                                          MatrixMemory,
                                          MatrixSize>,
                        Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                                      InputLayout,
                                      Kokkos::Cuda,
                                      InputMemory >,
                        Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                      OutputLayout,
                                      Kokkos::Cuda,
                                      OutputMemory >,
                        UpdateType
                      > {
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                             MatrixOrdinal,
                             Kokkos::Cuda,
                             MatrixMemory,
                             MatrixSize > MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                        InputLayout,
                        Kokkos::Cuda,
                        InputMemory > InputVectorType;
  typedef Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                        OutputLayout,
                        Kokkos::Cuda,
                        OutputMemory > OutputVectorType;
  typedef EnsembleMultiply<MatrixType,InputVectorType,OutputVectorType,UpdateType> type;
};

} // namespace details

} // namespace Stokhos

#endif /* #if defined( __CUDACC__) */

#endif /* #ifndef KOKKOS_CRSMATRIX_MP_VECTOR_CUDA_HPP */
