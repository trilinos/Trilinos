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
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"

#include "Kokkos_Blas1_UQ_PCE.hpp" // for some utilities

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsProductTensor.hpp"

namespace Stokhos {

//----------------------------------------------------------------------------
// Specialization of KokkosSparse::CrsMatrix for Sacado::UQ::PCE scalar type
//----------------------------------------------------------------------------

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef typename MatrixDevice::execution_space execution_space;

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename Kokkos::CijkType<matrix_values_type>::type tensor_type;
  typedef typename tensor_type::size_type size_type;
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

  Multiply( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y ,
            const input_scalar & a ,
            const output_scalar & b )
  : m_A_values( A.values )
  , m_A_graph( A.graph )
  , m_x( x )
  , m_y( y )
  , m_tensor( Kokkos::cijk(A.values) )
  , m_a( a )
  , m_b( b )
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
    if ( m_b == output_scalar(0) )
      for ( size_type j = 0 ; j < m_tensor.dimension() ; ++j )
        y[j] = 0 ;
    else
      for ( size_type j = 0 ; j < m_tensor.dimension() ; ++j )
        y[j] = m_b * y[j] ;
    //loop over cols of A
    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      const input_scalar * const x  = & m_x( 0 , m_A_graph.entries(iEntry) );
      const matrix_scalar * const A = & m_A_values( iEntry , 0 );

      BlockMultiply< tensor_type >::apply(
        m_tensor , A , x , y , m_a );
    }

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
  typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;
  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A_graph.row_map.extent(0)-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    const Kokkos::pair<size_type,size_type> work_range =
      details::compute_work_range<output_scalar>(
        device, m_tensor.dimension(), num_thread, thread_idx);

    // Prefer that y[ m_tensor.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    output_scalar * const y = & m_y(0,iBlockRow);

    // Leading dimension guaranteed contiguous for LayoutLeft
    if ( m_b == output_scalar(0) )
      for ( size_type j = work_range.first ; j < work_range.second ; ++j )
        y[j] = 0 ;
    else
      for ( size_type j = work_range.first ; j < work_range.second ; ++j )
        y[j] = m_b * y[j] ;

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

        // The number of nonzeros is always constrained to be a multiple of 8

        const size_type rem = iEntryEnd-iEntry;
        if (rem >= 8) {
          typedef TinyVec<tensor_scalar,8,tensor_type::use_intrinsics> ValTV2;
          typedef TinyVec<matrix_scalar,8,tensor_type::use_intrinsics> MatTV2;
          typedef TinyVec<output_scalar,8,tensor_type::use_intrinsics> VecTV2;
          const size_type *j = &m_tensor.coord(iEntry,0);
          const size_type *k = &m_tensor.coord(iEntry,1);
          ValTV2 c(&(m_tensor.value(iEntry)));

          for ( size_type col = 0; col < block_size; ++col ) {
            MatTV2 aj(sh_A[col], j), ak(sh_A[col], k);
            VecTV2 xj(sh_x[col], j), xk(sh_x[col], k);

            // vy += c * ( aj * xk + ak * xj)
            aj.times_equal(xk);
            aj.multiply_add(ak, xj);
            aj.times_equal(c);
            ytmp += aj.sum();
          }
        }

        y[iy] += m_a * ytmp ;
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
  typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;
  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A_graph.row_map.extent(0)-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    const Kokkos::pair<size_type,size_type> work_range =
      details::compute_work_range<output_scalar>(
        device, m_tensor.dimension(), num_thread, thread_idx);

    // Prefer that y[ m_tensor.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    output_scalar * const y = & m_y(0,iBlockRow);

    // Leading dimension guaranteed contiguous for LayoutLeft
    if ( m_b == output_scalar(0) )
      for ( size_type j = work_range.first ; j < work_range.second ; ++j )
        y[j] = 0 ;
    else
      for ( size_type j = work_range.first ; j < work_range.second ; ++j )
        y[j] = m_b * y[j] ;

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

        y[iy] += m_a * ytmp ;
      }

      // Add a team barrier to keep the thread team in-sync before going on
      // to the next block
      device.team_barrier();
    }

  }

#endif

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
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

    const size_t row_count = A.graph.row_map.extent(0) - 1 ;
    if (use_block_algorithm) {
#ifdef __MIC__
      const size_t team_size = 4;  // 4 hyperthreads for MIC
#else
      const size_t team_size = 2;  // 2 for everything else
#endif
      const size_t league_size = row_count;
      Kokkos::TeamPolicy< execution_space > config(league_size, team_size);
      Kokkos::parallel_for( config , Multiply(A,x,y,a,b) );
    }
    else {
      Kokkos::parallel_for( row_count , Multiply(A,x,y,a,b) );
    }
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>**,...>,
//   x and y are rank 2
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef typename MatrixDevice::execution_space execution_space;

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename Kokkos::CijkType<matrix_values_type>::type tensor_type;
  typedef typename tensor_type::size_type size_type;
  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
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

  Multiply( const matrix_type &        A ,
            const input_vector_type &  x ,
            const output_vector_type & y ,
            const input_scalar & a ,
            const output_scalar & b )
  : m_A_values( A.values )
  , m_A_graph( A.graph )
  , m_x( x )
  , m_y( y )
  , m_tensor( Kokkos::cijk(A.values) )
  , m_a( a )
  , m_b( b )
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

    const size_type num_col = m_y.extent(2);

    // Leading dimension guaranteed contiguous for LayoutLeft
    if ( m_b == output_scalar(0) )
      for (size_type col=0; col<num_col; ++col)
        for ( size_type j = 0 ; j < m_tensor.dimension() ; ++j )
          m_y(j, iBlockRow, col) = 0 ;
    else
       for (size_type col=0; col<num_col; ++col)
        for ( size_type j = 0 ; j < m_tensor.dimension() ; ++j )
          m_y(j, iBlockRow, col) = m_b * m_y(j, iBlockRow, col) ;

    // Put the x-column loop inside the A-column loop to reuse entries in A.
    // This way all of the entries for that particular column of A should stay
    // in L1 cache for all of the columns of x.

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      const matrix_scalar * const A = &m_A_values( iEntry, 0 );
      const size_type iBlockCol = m_A_graph.entries(iEntry);

      for (size_type col=0; col<num_col; ++col) {
        output_scalar * const y =      &m_y( 0, iBlockRow, col );
        const input_scalar * const x = &m_x( 0, iBlockCol, col );
        BlockMultiply< tensor_type >::apply( m_tensor , A , x , y , m_a );
      }

    }

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
  typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A_graph.row_map.extent(0)-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    const Kokkos::pair<size_type,size_type> work_range =
      details::compute_work_range<output_scalar>(
        device, m_tensor.dimension(), num_thread, thread_idx);

    const size_type num_col = m_y.extent(2);

    // Leading dimension guaranteed contiguous for LayoutLeft
    if ( m_b == output_scalar(0) )
      for (size_type col=0; col<num_col; ++col)
        for ( size_type j = work_range.first ; j < work_range.second ; ++j )
          m_y(j, iBlockRow, col) = 0 ;
    else
      for (size_type col=0; col<num_col; ++col)
        for ( size_type j = work_range.first ; j < work_range.second ; ++j )
           m_y(j, iBlockRow, col) = m_b *  m_y(j, iBlockRow, col) ;

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

      // Loop over columns of x, y
      for (size_type vec_col=0; vec_col<num_col; ++vec_col) {

        output_scalar * const y = & m_y( 0 , iBlockRow , vec_col );

        for ( size_type col = 0; col < block_size; ++col ) {
          const size_type iBlockColumn = m_A_graph.entries( iBlockEntry + col );
          sh_x[col] = & m_x( 0 , iBlockColumn  , vec_col );
          sh_A[col] = & m_A_values( iBlockEntry + col , 0);
        }

        for ( size_type iy = work_range.first ; iy < work_range.second ; ++iy ){

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

          // The number of nonzeros is always constrained to be a multiple of 8

          const size_type rem = iEntryEnd-iEntry;
          if (rem >= 8) {
            typedef TinyVec<tensor_scalar,8,tensor_type::use_intrinsics> ValTV2;
            typedef TinyVec<matrix_scalar,8,tensor_type::use_intrinsics> MatTV2;
            typedef TinyVec<output_scalar,8,tensor_type::use_intrinsics> VecTV2;
            const size_type *j = &m_tensor.coord(iEntry,0);
            const size_type *k = &m_tensor.coord(iEntry,1);
            ValTV2 c(&(m_tensor.value(iEntry)));

            for ( size_type col = 0; col < block_size; ++col ) {
              MatTV2 aj(sh_A[col], j), ak(sh_A[col], k);
              VecTV2 xj(sh_x[col], j), xk(sh_x[col], k);

              // vy += c * ( aj * xk + ak * xj)
              aj.times_equal(xk);
              aj.multiply_add(ak, xj);
              aj.times_equal(c);
              ytmp += aj.sum();
            }
          }

          y[iy] += m_a * ytmp ;
        }

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
  typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A_graph.row_map.extent(0)-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    const Kokkos::pair<size_type,size_type> work_range =
      details::compute_work_range<output_scalar>(
        device, m_tensor.dimension(), num_thread, thread_idx);

    const size_type num_col = m_y.extent(2);

    // Leading dimension guaranteed contiguous for LayoutLeft
    if ( m_b == output_scalar(0) )
      for (size_type col=0; col<num_col; ++col)
        for ( size_type j = work_range.first ; j < work_range.second ; ++j )
          m_y(j, iBlockRow, col) = 0 ;
    else
      for (size_type col=0; col<num_col; ++col)
        for ( size_type j = work_range.first ; j < work_range.second ; ++j )
           m_y(j, iBlockRow, col) = m_b *  m_y(j, iBlockRow, col) ;

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

      // Loop over columns of x, y
      for (size_type vec_col=0; vec_col<num_col; ++vec_col) {

        output_scalar * const y = & m_y( 0 , iBlockRow , vec_col );

        for ( size_type col = 0; col < block_size; ++col ) {
          const size_type iBlockColumn = m_A_graph.entries( iBlockEntry + col );
          sh_x[col] = & m_x( 0 , iBlockColumn , vec_col );
          sh_A[col] = & m_A_values( iBlockEntry + col , 0 );
        }

        for ( size_type iy = work_range.first ; iy < work_range.second ; ++iy ){

          const size_type nEntry = m_tensor.num_entry(iy);
          const size_type iEntryBeg = m_tensor.entry_begin(iy);
          const size_type iEntryEnd = iEntryBeg + nEntry;
                size_type iEntry    = iEntryBeg;

          output_scalar ytmp = 0 ;

          // Do entries with a blocked loop of size blocksize
          if (tensor_type::vectorsize > 1 && nEntry >= tensor_type::vectorsize){
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

          y[iy] += m_a * ytmp ;
        }

      }

      // Add a team barrier to keep the thread team in-sync before going on
      // to the next block
      device.team_barrier();
    }

  }

#endif

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
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

    const size_t row_count = A.graph.row_map.extent(0) - 1 ;
    if (use_block_algorithm) {
#ifdef __MIC__
      const size_t team_size = 4;  // 4 hyperthreads for MIC
#else
      const size_t team_size = 2;  // 2 for everything else
#endif
      const size_t league_size = row_count;
      Kokkos::TeamPolicy< execution_space > config(league_size, team_size);
      Kokkos::parallel_for( config , Multiply(A,x,y,a,b) );
    }
    else {
      Kokkos::parallel_for( row_count , Multiply(A,x,y,a,b) );
    }
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef typename matrix_type::const_type const_matrix_type;

  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;

  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const_matrix_type cA = A;
    Multiply< const_matrix_type, input_vector_type, output_vector_type >::apply(
      cA, x, y, a, b);
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>**,...>,
//   x and y are rank 2
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef typename matrix_type::const_type const_matrix_type;

  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;

  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const_matrix_type cA = A;
    Multiply< const_matrix_type, input_vector_type, output_vector_type >::apply(
      cA, x, y, a, b);
  }
};

template <typename MatrixType, typename InputViewType, typename OutputViewType>
class MeanMultiply {};

// Kernel implementing y = A * x where PCE size of A is 1
//   A == KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<...>,...>, with A.values.sacado_size() == 1
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class MeanMultiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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
  struct BlockKernel {
    typedef typename MatrixDevice::execution_space execution_space;
    typedef typename Kokkos::FlatArrayType<matrix_values_type>::type matrix_array_type;
    typedef typename input_vector_type::array_type input_array_type;
    typedef typename output_vector_type::array_type output_array_type;

    const matrix_array_type   m_A_values ;
    const matrix_graph_type   m_A_graph ;
    const output_array_type   v_y ;
    const input_array_type    v_x ;
    const input_scalar        m_a ;
    const output_scalar       m_b ;
    const size_type           dim ;
    const size_type           numBlocks ;
    const size_type           rem ;
    const size_type           dim_block ;

    BlockKernel( const matrix_type &        A ,
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
      , dim( dimension_scalar(x) )
      , numBlocks( dim / BlockSize )
      , rem( dim % BlockSize )
      , dim_block( numBlocks*BlockSize )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type iBlockRow ) const
    {
#if defined(__INTEL_COMPILER)&& ! defined(__CUDA_ARCH__)
      output_scalar s[BlockSize] __attribute__((aligned(64))) = {};
#else
      output_scalar s[BlockSize] = {};
#endif

      const size_type iEntryBegin = m_A_graph.row_map[ iBlockRow ];
      const size_type iEntryEnd   = m_A_graph.row_map[ iBlockRow + 1 ];
      size_type pce_block = 0;
      for (; pce_block < dim_block; pce_block+=BlockSize) {
        output_scalar * const y = &v_y(pce_block, iBlockRow);
        if (m_b == output_scalar(0))
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
          for (size_type k = 0; k < BlockSize; ++k)
            s[k] = 0.0;
        else
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
          for (size_type k = 0; k < BlockSize; ++k)
            s[k] = m_b*y[k];
        for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
          const matrix_scalar aA = m_a*m_A_values(iEntry);
          const size_type col = m_A_graph.entries(iEntry);
          const input_scalar * const x = &v_x(pce_block, col);
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
          for (size_type k = 0; k < BlockSize; ++k)
            s[k] += aA*x[k];
        }
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
        for (size_type k = 0; k < BlockSize; ++k) {
          y[k] = s[k];
        }
      }

      // Remaining coeffs
      if (rem > 0) {
        output_scalar * const y = &v_y(pce_block, iBlockRow);
        if (m_b == output_scalar(0))
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
          for (size_type k = 0; k < rem; ++k)
            s[k] = 0.0;
        else
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
          for (size_type k = 0; k < rem; ++k)
            s[k] = m_b*y[k];
        for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
          const matrix_scalar aA = m_a*m_A_values(iEntry);
          const size_type col = m_A_graph.entries(iEntry);
          const input_scalar * const x = &v_x(pce_block, col);
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
          for (size_type k = 0; k < rem; ++k)
            s[k] += aA*x[k];
        }
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
        for (size_type k = 0; k < rem; ++k) {
          y[k] = s[k];
        }
      }

    }

  };

  struct Kernel {
    typedef typename MatrixDevice::execution_space execution_space;
    typedef typename Kokkos::FlatArrayType<matrix_values_type>::type matrix_array_type;
    typedef typename input_vector_type::array_type input_array_type;
    typedef typename output_vector_type::array_type output_array_type;

    const matrix_array_type   m_A_values ;
    const matrix_graph_type   m_A_graph ;
    const output_array_type   v_y ;
    const input_array_type    v_x ;
    const input_scalar        m_a ;
    const output_scalar       m_b ;
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
      , dim( dimension_scalar(x) )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type iBlockRow ) const
    {
      const size_type iEntryBegin = m_A_graph.row_map[ iBlockRow ];
      const size_type iEntryEnd   = m_A_graph.row_map[ iBlockRow + 1 ];
      output_scalar * const y = &v_y(0, iBlockRow);
      if (m_b == output_scalar(0))
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
        for (size_type k = 0; k < dim; ++k)
          y[k] = 0.0;
      else
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
        for (size_type k = 0; k < dim; ++k)
          y[k] = m_b*y[k];
      for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
        const matrix_scalar aA = m_a*m_A_values(iEntry);
        const size_type col = m_A_graph.entries(iEntry);
        const input_scalar * const x = &v_x(0, col);
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma ivdep
//#pragma vector aligned
#endif
        for (size_type k = 0; k < dim; ++k)
          y[k] += aA*x[k];
      }
    }

  };

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const size_t row_count = A.graph.row_map.extent(0) - 1 ;
    const size_type dim = Kokkos::dimension_scalar(x);

    // Choose block size appropriately for PCE dimension
#if defined (__MIC__)
    if (dim >= 128)
      Kokkos::parallel_for( row_count , Kernel(A,x,y,a,b) );
    else if (dim >= 64)
      Kokkos::parallel_for( row_count , BlockKernel<64>(A,x,y,a,b) );
    else if (dim >= 32)
      Kokkos::parallel_for( row_count , BlockKernel<32>(A,x,y,a,b) );
    else if (dim >= 16)
      Kokkos::parallel_for( row_count , BlockKernel<16>(A,x,y,a,b) );
    else if (dim >= 8)
      Kokkos::parallel_for( row_count , BlockKernel<8>(A,x,y,a,b) );
    else
      Kokkos::parallel_for( row_count , BlockKernel<4>(A,x,y,a,b) );
#else
    if (dim >= 32)
      Kokkos::parallel_for( row_count , BlockKernel<32>(A,x,y,a,b) );
    else if (dim >= 16)
      Kokkos::parallel_for( row_count , BlockKernel<16>(A,x,y,a,b) );
    else if (dim >= 8)
      Kokkos::parallel_for( row_count , BlockKernel<8>(A,x,y,a,b) );
    else
      Kokkos::parallel_for( row_count , BlockKernel<4>(A,x,y,a,b) );
#endif
  }
};

// Kernel implementing y = A * x where A has PCE size = 1
//   A == KokkosSparse::CrsMatrix< Sacado::UQ::PCE<...>,...>,
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>**,...>,
//   x and y are rank 2
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class MeanMultiply< KokkosSparse::CrsMatrix< const Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename MatrixValue::ordinal_type size_type;
  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    typedef Kokkos::View< const InputVectorValue*,
      InputP... > input_vector_1d_type;
    typedef Kokkos::View< OutputVectorValue*,
      OutputP... > output_vector_1d_type;
    typedef MeanMultiply<matrix_type, input_vector_1d_type,
      output_vector_1d_type> MeanMultiply1D;
    const size_type num_col = x.extent(1);
    for (size_type i=0; i<num_col; ++i) {
      input_vector_1d_type x_col =
        Kokkos::subview(x, Kokkos::ALL(), i);
      output_vector_1d_type y_col =
        Kokkos::subview(y, Kokkos::ALL(), i);
      MeanMultiply1D::apply( A, x_col, y_col, a, b );
    }
  }
};

// Kernel implementing y = A * x where PCE size of A is 1
//   A == KokkosSparse::CrsMatrix< Sacado::UQ::PCE<...>,...>, with A.values.sacado_size() == 1
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>*,...>,
//   x and y are rank 1
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class MeanMultiply< KokkosSparse::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef typename matrix_type::const_type const_matrix_type;

  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;

  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const_matrix_type cA = A;
    MeanMultiply< const_matrix_type, input_vector_type, output_vector_type >::apply(
      cA, x, y, a, b);
  }
};

// Kernel implementing y = A * x where PCE size of A is 1
//   A == KokkosSparse::CrsMatrix< Sacado::UQ::PCE<...>,...>, with A.values.sacado_size() == 1
//   x, y == Kokkos::View< Sacado::UQ::PCE<...>**,...>,
//   x and y are rank 2
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class MeanMultiply< KokkosSparse::CrsMatrix< Sacado::UQ::PCE<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
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

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize> matrix_type;
  typedef typename matrix_type::const_type const_matrix_type;

  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;

  typedef typename InputVectorValue::value_type input_scalar;
  typedef typename OutputVectorValue::value_type output_scalar;

  static void apply( const matrix_type & A ,
                     const input_vector_type & x ,
                     const output_vector_type & y ,
                     const input_scalar & a = input_scalar(1) ,
                     const output_scalar & b = output_scalar(0) )
  {
    const_matrix_type cA = A;
    MeanMultiply< const_matrix_type, input_vector_type, output_vector_type >::apply(
      cA, x, y, a, b);
  }
};

} // namespace Stokhos

namespace KokkosSparse {

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_uq_pce< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_ONE)
{
  typedef Kokkos::View< OutputType, OutputP... > OutputVectorType;
  typedef Kokkos::View< InputType, InputP... > InputVectorType;
  typedef Stokhos::Multiply<MatrixType, typename InputVectorType::const_type,
                            OutputVectorType> multiply_type;
  typedef Stokhos::MeanMultiply<MatrixType, typename InputVectorType::const_type,
                                OutputVectorType> mean_multiply_type;

  if(mode[0]!='N') {
    Kokkos::Impl::raise_error(
      "Stokhos spmv not implemented for transposed or conjugated matrix-vector multiplies");
  }

  if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
    Kokkos::Impl::raise_error(
      "Stokhos spmv not implemented for non-constant a or b");
  }
  if (dimension_scalar(A.values) == 1 && dimension_scalar(x) != 1) {
    mean_multiply_type::apply( A, x, y,
                               Sacado::Value<AlphaType>::eval(a),
                               Sacado::Value<BetaType>::eval(b) );
  }
  else
    multiply_type::apply( A, x, y,
                          Sacado::Value<AlphaType>::eval(a),
                          Sacado::Value<BetaType>::eval(b) );
}

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_uq_pce< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  KokkosKernels::Experimental::Controls,
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_ONE)
{
  spmv(mode, a, A, x, b, y, RANK_ONE());
}

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_uq_pce< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_TWO)
{
  if(mode[0]!='N') {
    Kokkos::Impl::raise_error(
      "Stokhos spmv not implemented for transposed or conjugated matrix-vector multiplies");
  }
  if (y.extent(1) == 1) {
    auto y_1D = subview(y, Kokkos::ALL(), 0);
    auto x_1D = subview(x, Kokkos::ALL(), 0);
    spmv(mode, a, A, x_1D, b, y_1D, RANK_ONE());
  }
  else {
    typedef Kokkos::View< OutputType, OutputP... > OutputVectorType;
    typedef Kokkos::View< InputType, InputP... > InputVectorType;

    typedef Stokhos::Multiply<MatrixType, typename InputVectorType::const_type,
                              OutputVectorType> multiply_type;
    typedef Stokhos::MeanMultiply<MatrixType, typename InputVectorType::const_type,
                                  OutputVectorType> mean_multiply_type;

    if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
      Kokkos::Impl::raise_error(
        "Stokhos spmv not implemented for non-constant a or b");
    }

    typename InputVectorType::const_type x_const = x;

    if (dimension_scalar(A.values) == 1 && dimension_scalar(x) != 1) {
      mean_multiply_type::apply( A, x_const, y,
                                 Sacado::Value<AlphaType>::eval(a),
                                 Sacado::Value<BetaType>::eval(b));
     }
    else
      multiply_type::apply( A, x_const, y,
                            Sacado::Value<AlphaType>::eval(a),
                            Sacado::Value<BetaType>::eval(b));
  }
}

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_uq_pce< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  KokkosKernels::Experimental::Controls,
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_TWO)
{
  spmv(mode, a, A, x, b, y, RANK_TWO());
}

}

#endif /* #ifndef KOKKOS_CRSMATRIX_UQ_PCE_HPP */
