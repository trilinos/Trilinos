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

#ifndef KOKKOS_CRSMATRIX_UQ_PCE_HPP
#define KOKKOS_CRSMATRIX_UQ_PCE_HPP

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_MV_UQ_PCE.hpp" // for some utilities

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsProductTensor.hpp"

namespace Stokhos {

//----------------------------------------------------------------------------
// Specialization of Kokkos CrsMatrix math functions
//----------------------------------------------------------------------------

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputMemory,
          typename OutputStorage,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::UQ::PCE<InputStorage>*,
                              Kokkos::LayoutLeft,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::UQ::PCE<OutputStorage>*,
                              Kokkos::LayoutLeft,
                              Device,
                              OutputMemory >
                >
{
public:
  typedef Sacado::UQ::PCE<MatrixStorage> MatrixValue;
  typedef Sacado::UQ::PCE<InputStorage> InputVectorValue;
  typedef Sacado::UQ::PCE<OutputStorage> OutputVectorValue;

  typedef Device device_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_values_type::cijk_type tensor_type;
  typedef typename tensor_type::size_type size_type;
  typedef Kokkos::View< InputVectorValue*,
                        Kokkos::LayoutLeft,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        Kokkos::LayoutLeft,
                        Device,
                        OutputMemory > output_vector_type;

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

  Multiply( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y )
  : m_A_values( A.values )
  , m_A_graph( A.graph )
  , m_x( x )
  , m_y( y )
  , m_tensor( A.values.cijk() )
  {}

public:

  //
  // Non-team functor interface -- no threads within PCE multiply
  //
  // Note:  Rember that matrix currently is always LayoutRight!
  //
  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iBlockRow ) const
  {
    // Prefer that y[ m_tensor.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    output_scalar * const y = & m_y(0,iBlockRow);

    const size_type iEntryBegin = m_A_graph.row_map[ iBlockRow ];
    const size_type iEntryEnd   = m_A_graph.row_map[ iBlockRow + 1 ];

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = 0 ; j < m_tensor.dimension() ; ++j ) { y[j] = 0 ; }

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      const input_scalar * const x = & m_x( 0 , m_A_graph.entries(iEntry) );
      const matrix_scalar * const a = & m_A_values( iEntry , 0 );

      BlockMultiply< tensor_type >::apply( m_tensor , a , x , y );
    }

  }

  /*
   * Compute work range = (begin, end) such that adjacent threads write to
   * separate cache lines
   */
  KOKKOS_INLINE_FUNCTION
  std::pair< size_type , size_type >
  compute_work_range( const size_type work_count ,
                      const size_type thread_count ,
                      const size_type thread_rank ) const
  {
    enum { work_align = 64 / sizeof(output_scalar) };
    enum { work_shift = Kokkos::Impl::power_of_two< work_align >::value };
    enum { work_mask  = work_align - 1 };

    const size_type work_per_thread =
      ( ( ( ( work_count + work_mask ) >> work_shift ) + thread_count - 1 ) /
        thread_count ) << work_shift ;

    const size_type work_begin =
      std::min( thread_rank * work_per_thread , work_count );
    const size_type work_end   =
      std::min( work_begin + work_per_thread , work_count );

    return std::make_pair( work_begin , work_end );
  }

#if defined(__MIC__)

  //
  // Team functor interface with threading within PCE multiply
  //
  // Note:  Rember that matrix currently is always LayoutRight!
  //
  // This is a MIC-specific version of that processes multiple FEM columns
  // at a time to reduce tensor reads
  //
  KOKKOS_INLINE_FUNCTION
  void operator()( device_type device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A_graph.row_map.dimension_0()-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    std::pair<size_type,size_type> work_range =
      compute_work_range(m_tensor.dimension(), num_thread, thread_idx);

    // Prefer that y[ m_tensor.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    output_scalar * const y = & m_y(0,iBlockRow);

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = work_range.first ; j < work_range.second ; ++j )
      y[j] = 0 ;

    const size_type iBlockEntryBeg = m_A_graph.row_map[ iBlockRow ];
    const size_type iBlockEntryEnd = m_A_graph.row_map[ iBlockRow + 1 ];
    const size_type BlockSize = 9;
    const size_type numBlock =
      (iBlockEntryEnd-iBlockEntryBeg+BlockSize-1) / BlockSize;

    const matrix_scalar* sh_A[BlockSize];
    const input_scalar* sh_x[BlockSize];

    size_type iBlockEntry = iBlockEntryBeg;
    for (size_type block = 0; block<numBlock; ++block, iBlockEntry+=BlockSize) {
      const size_type block_size =
        block == numBlock-1 ? iBlockEntryEnd-iBlockEntry : BlockSize;

      for ( size_type col = 0; col < block_size; ++col ) {
        const size_type iBlockColumn = m_A_graph.entries( iBlockEntry + col );
        sh_x[col] = & m_x( 0 , iBlockColumn );
        sh_A[col] = & m_A_values( iBlockEntry + col , 0);
      }

      for ( size_type iy = work_range.first ; iy < work_range.second ; ++iy ) {

        const size_type nEntry = m_tensor.num_entry(iy);
        const size_type iEntryBeg = m_tensor.entry_begin(iy);
        const size_type iEntryEnd = iEntryBeg + nEntry;
              size_type iEntry    = iEntryBeg;

        output_scalar ytmp = 0 ;

        // Do entries with a blocked loop of size blocksize
        const size_type nBlock = nEntry / tensor_type::vectorsize;
        const size_type nEntryB = nBlock * tensor_type::vectorsize;
        const size_type iEnd = iEntryBeg + nEntryB;

        typedef TinyVec<tensor_scalar,tensor_type::vectorsize,tensor_type::use_intrinsics> ValTV;
        typedef TinyVec<matrix_scalar,tensor_type::vectorsize,tensor_type::use_intrinsics> MatTV;
        typedef TinyVec<output_scalar,tensor_type::vectorsize,tensor_type::use_intrinsics> VecTV;
        VecTV vy;
        vy.zero();
        for (size_type block=0; block<nBlock; ++block, iEntry+=tensor_type::vectorsize) {
          const size_type *j = &m_tensor.coord(iEntry,0);
          const size_type *k = &m_tensor.coord(iEntry,1);
          ValTV c(&(m_tensor.value(iEntry)));

          for ( size_type col = 0; col < block_size; ++col ) {
            MatTV aj(sh_A[col], j), ak(sh_A[col], k);
            VecTV xj(sh_x[col], j), xk(sh_x[col], k);

            // vy += c * ( aj * xk + ak * xj)
            aj.times_equal(xk);
            aj.multiply_add(ak, xj);
            vy.multiply_add(c, aj);
          }
        }
        ytmp += vy.sum();

        const size_type rem = iEntryEnd-iEntry;
        if (rem > 8) {
          typedef TinyVec<tensor_scalar,tensor_type::vectorsize,true,true> ValTV2;
          typedef TinyVec<matrix_scalar,tensor_type::vectorsize,true,true> MatTV2;
          typedef TinyVec<output_scalar,tensor_type::vectorsize,true,true> VecTV2;
          const size_type *j = &m_tensor.coord(iEntry,0);
          const size_type *k = &m_tensor.coord(iEntry,1);
          ValTV2 c(&(m_tensor.value(iEntry)), rem);

          for ( size_type col = 0; col < block_size; ++col ) {
            MatTV2 aj(sh_A[col], j, rem), ak(sh_A[col], k, rem);
            VecTV2 xj(sh_x[col], j, rem), xk(sh_x[col], k, rem);

            // vy += c * ( aj * xk + ak * xj)
            aj.times_equal(xk);
            aj.multiply_add(ak, xj);
            aj.times_equal(c);
            ytmp += aj.sum();
            iEntry += rem;
          }
        }

        else if (rem > 0) {
          typedef TinyVec<tensor_scalar,8,true,true> ValTV2;
          typedef TinyVec<matrix_scalar,8,true,true> MatTV2;
          typedef TinyVec<output_scalar,8,true,true> VecTV2;
          const size_type *j = &m_tensor.coord(iEntry,0);
          const size_type *k = &m_tensor.coord(iEntry,1);
          ValTV2 c(&(m_tensor.value(iEntry)), rem);

          for ( size_type col = 0; col < block_size; ++col ) {
            MatTV2 aj(sh_A[col], j, rem), ak(sh_A[col], k, rem);
            VecTV2 xj(sh_x[col], j, rem), xk(sh_x[col], k, rem);

            // vy += c * ( aj * xk + ak * xj)
            aj.times_equal(xk);
            aj.multiply_add(ak, xj);
            aj.times_equal(c);
            ytmp += aj.sum();
            iEntry += rem;
          }
        }

        y[iy] += ytmp ;
      }

      // Add a team barrier to keep the thread team in-sync before going on
      // to the next block
      device.team_barrier();
    }

  }

#else

  //
  // Team functor interface with threading within PCE multiply
  //
  // Note:  Rember that matrix currently is always LayoutRight!
  //
  // This is a general, hand-vectorized version that processes multiple FEM
  // columns at a time to reduce tensor reads.  Note that auto-vectorization
  // doesn't work here because of the inner-loop over FEM columns.
  //
  KOKKOS_INLINE_FUNCTION
  void operator()( device_type device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A_graph.row_map.dimension_0()-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    std::pair<size_type,size_type> work_range =
      compute_work_range(m_tensor.dimension(), num_thread, thread_idx);

    // Prefer that y[ m_tensor.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    output_scalar * const y = & m_y(0,iBlockRow);

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = work_range.first ; j < work_range.second ; ++j )
      y[j] = 0 ;

    const size_type iBlockEntryBeg = m_A_graph.row_map[ iBlockRow ];
    const size_type iBlockEntryEnd = m_A_graph.row_map[ iBlockRow + 1 ];
    const size_type BlockSize = 14;
    const size_type numBlock =
      (iBlockEntryEnd-iBlockEntryBeg+BlockSize-1) / BlockSize;

    const matrix_scalar* sh_A[BlockSize];
    const input_scalar* sh_x[BlockSize];

    size_type iBlockEntry = iBlockEntryBeg;
    for (size_type block = 0; block<numBlock; ++block, iBlockEntry+=BlockSize) {
      const size_type block_size =
        block == numBlock-1 ? iBlockEntryEnd-iBlockEntry : BlockSize;

      for ( size_type col = 0; col < block_size; ++col ) {
        const size_type iBlockColumn = m_A_graph.entries( iBlockEntry + col );
        sh_x[col] = & m_x( 0 , iBlockColumn );
        sh_A[col] = & m_A_values( iBlockEntry + col , 0 );
      }

      for ( size_type iy = work_range.first ; iy < work_range.second ; ++iy ) {

        const size_type nEntry = m_tensor.num_entry(iy);
        const size_type iEntryBeg = m_tensor.entry_begin(iy);
        const size_type iEntryEnd = iEntryBeg + nEntry;
              size_type iEntry    = iEntryBeg;

        output_scalar ytmp = 0 ;

        // Do entries with a blocked loop of size blocksize
        if (tensor_type::vectorsize > 1 && nEntry >= tensor_type::vectorsize) {
          const size_type nBlock = nEntry / tensor_type::vectorsize;
          const size_type nEntryB = nBlock * tensor_type::vectorsize;
          const size_type iEnd = iEntryBeg + nEntryB;

          typedef TinyVec<tensor_scalar,tensor_type::vectorsize,tensor_type::use_intrinsics> ValTV;
          typedef TinyVec<matrix_scalar,tensor_type::vectorsize,tensor_type::use_intrinsics> MatTV;
          typedef TinyVec<output_scalar,tensor_type::vectorsize,tensor_type::use_intrinsics> VecTV;
          VecTV vy;
          vy.zero();
          for (; iEntry<iEnd; iEntry+=tensor_type::vectorsize) {
            const size_type *j = &m_tensor.coord(iEntry,0);
            const size_type *k = &m_tensor.coord(iEntry,1);
            ValTV c(&(m_tensor.value(iEntry)));

            for ( size_type col = 0; col < block_size; ++col ) {
              MatTV aj(sh_A[col], j), ak(sh_A[col], k);
              VecTV xj(sh_x[col], j), xk(sh_x[col], k);

              // vy += c * ( aj * xk + ak * xj)
              aj.times_equal(xk);
              aj.multiply_add(ak, xj);
              vy.multiply_add(c, aj);
            }
          }
          ytmp += vy.sum();
        }

        // Do remaining entries with a scalar loop
        for ( ; iEntry<iEntryEnd; ++iEntry) {
          const size_type j = m_tensor.coord(iEntry,0);
          const size_type k = m_tensor.coord(iEntry,1);
          tensor_scalar cijk = m_tensor.value(iEntry);

          for ( size_type col = 0; col < block_size; ++col ) {
            ytmp += cijk * ( sh_A[col][j] * sh_x[col][k] +
                             sh_A[col][k] * sh_x[col][j] );
          }

        }

        y[iy] += ytmp ;
      }

      // Add a team barrier to keep the thread team in-sync before going on
      // to the next block
      device.team_barrier();
    }

  }

#endif

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y )
  {
    // Generally the block algorithm seems to perform better on the MIC,
    // as long as the stochastic size isn't too big, but doesn't perform
    // any better on the CPU (probably because the CPU has a fat L3 cache
    // to store the sparse 3 tensor).
#ifdef __MIC__
    const bool use_block_algorithm = true;
#else
    const bool use_block_algorithm = false;
#endif

    const size_t row_count = A.graph.row_map.dimension_0() - 1 ;
    if (use_block_algorithm) {
      typedef typename matrix_type::device_type device_type;
#ifdef __MIC__
      const size_t team_size = 4;  // 4 hyperthreads for MIC
#else
      const size_t team_size = 2;  // 2 for everything else
#endif
      const size_t league_size = row_count;
      Kokkos::ParallelWorkRequest config(league_size, team_size);
      Kokkos::parallel_for( config , Multiply(A,x,y) );
    }
    else {
      Kokkos::parallel_for( row_count , Multiply(A,x,y) );
    }
  }
};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>**,...>,
//   x and y are rank 2
//
// Note:  Unlike the rank-1 version, this version has not been
// optimized, and doesn't even include the block-column implementation
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputMemory,
          typename OutputStorage,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::UQ::PCE<InputStorage>**,
                              Kokkos::LayoutLeft,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::UQ::PCE<OutputStorage>**,
                              Kokkos::LayoutLeft,
                              Device,
                              OutputMemory >
                >
{
public:
  typedef Sacado::UQ::PCE<MatrixStorage> MatrixValue;
  typedef Sacado::UQ::PCE<InputStorage> InputVectorValue;
  typedef Sacado::UQ::PCE<OutputStorage> OutputVectorValue;

  typedef Device device_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_values_type::cijk_type tensor_type;
  typedef typename tensor_type::size_type size_type;
  typedef Kokkos::View< InputVectorValue**,
                        Kokkos::LayoutLeft,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        Kokkos::LayoutLeft,
                        Device,
                        OutputMemory > output_vector_type;

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

  Multiply( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y )
  : m_A_values( A.values )
  , m_A_graph( A.graph )
  , m_x( x )
  , m_y( y )
  , m_tensor( A.values.cijk() )
  {}

public:

  //
  // Non-team functor interface -- no threads within PCE multiply
  //
  // Note:  Rember that matrix currently is always LayoutRight!
  //
  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iBlockRow ) const
  {
    const size_type iEntryBegin = m_A_graph.row_map[ iBlockRow ];
    const size_type iEntryEnd   = m_A_graph.row_map[ iBlockRow + 1 ];

    const size_type num_col = m_y.dimension_2();

    // Leading dimension guaranteed contiguous for LayoutLeft
    for (size_type col=0; col<num_col; ++col)
      for ( size_type j = 0 ; j < m_tensor.dimension() ; ++j )
        m_y(j, iBlockRow, col) = 0 ;

    // Put the x-column loop inside the A-column loop to reuse entries in A.
    // This way all of the entries for that particular column of A should stay
    // in L1 cache for all of the columns of x.

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      const matrix_scalar * const a = &m_A_values( iEntry, 0 );
      const size_type iBlockCol = m_A_graph.entries(iEntry);

      for (size_type col=0; col<num_col; ++col) {
        output_scalar * const y =      &m_y( 0, iBlockRow, col );
        const input_scalar * const x = &m_x( 0, iBlockCol, col );
        BlockMultiply< tensor_type >::apply( m_tensor , a , x , y );
      }

    }

  }

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y )
  {
    const size_t row_count = A.graph.row_map.dimension_0() - 1 ;
    Kokkos::parallel_for( row_count , Multiply(A,x,y) );
  }
};

} // namespace Stokhos

namespace Kokkos {

// Overload of Kokkos::MV_Multiply for Sacado::UQ::PCE scalar types
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
  const Kokkos::View< Sacado::UQ::PCE< OutputStorage>*,
                OutputLayout,
                Device,
                OutputMemory >& y,
  const Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::UQ::PCE<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef Kokkos::View< Sacado::UQ::PCE< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::UQ::PCE<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;
  typedef Stokhos::Multiply<MatrixType,InputVectorType,
    OutputVectorType> multiply_type;
  multiply_type::apply( A, x, y );
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
  const Kokkos::View< Sacado::UQ::PCE< OutputStorage>*,
                OutputLayout,
                Device,
                OutputMemory >& y,
  const Sacado::UQ::PCE<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::UQ::PCE<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef typename InputStorage::value_type value_type;
  if (Impl::is_pce_constant(a) && a.fastAccessCoeff(0) == value_type(1)) {
    MV_Multiply(y, A, x);
  }
  else {
    Impl::raise_sacado_error(
      "MV_Multiply not implemented for non-constant a != 1");
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
  const Kokkos::View< Sacado::UQ::PCE< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Sacado::UQ::PCE<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::UQ::PCE<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::UQ::PCE< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::UQ::PCE<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(y_1D, a, A, x_1D);
  }
  else {
    typedef typename InputStorage::value_type value_type;
    if (Impl::is_pce_constant(a) && a.fastAccessCoeff(0) == value_type(1)) {
      typedef Kokkos::View< Sacado::UQ::PCE< OutputStorage>**,
        OutputLayout, Device, OutputMemory > OutputVectorType;
      typedef Kokkos::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
        MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
      typedef Kokkos::View< Sacado::UQ::PCE<InputStorage>**,
        InputLayout, Device, InputMemory > InputVectorType;
      typedef Stokhos::Multiply<MatrixType,InputVectorType,
        OutputVectorType> multiply_type;
      multiply_type::apply( A, x, y );
      //Stokhos::multiply(A, x, y, Tag());
    }
    else {
      Impl::raise_sacado_error(
        "MV_Multiply not implemented for non-constant a != 1");
    }
  }
}

}

#endif /* #ifndef KOKKOS_CRSMATRIX_UQ_PCE_HPP */
