/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

//#include "KokkosKernels_ETIHelperMacros.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <stdexcept>
#include "KokkosSparse_BlockCrsMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

namespace Test{ // anonymous

  using std::cerr;
  using std::endl;

  // Create a test sparse matrix A.
  //
  // Identify the matrix to create by number (whichMatrix).  The
  // following lists the valid options for whichMatrix:
  //
  // 0: A square 8 x 8 sparse CrsMatrix with implicit block structure
  // 1: A square 4 x 4 sparse BlockCrsMatrix
  //
  // \param ptr [out] Array of row offsets, of length numRows+1.
  // \param ind [out] Array of column indices, of length nnz (CrsMatrix) 
  //        or numBlocks (BlockCrsMatrix).
  // \param val [out] Array of entries (values), of length nnz.
  // \param numRows [out] The number of rows in the matrix.
  // \param numCols [out] The number of columns in the matrix.
  // \param nnz [out] The number of stored entries in the matrix.
  // \param whichMatrix [in] The index of the matrix to create.
  template<typename sparseMat_t>
  void
  makeSparseMatrix (
      typename sparseMat_t::StaticCrsGraphType::row_map_type::non_const_type & ptr,
      typename sparseMat_t::StaticCrsGraphType::entries_type::non_const_type   & ind,
      typename sparseMat_t::values_type::non_const_type & val,
      typename sparseMat_t::ordinal_type &numRows,
      typename sparseMat_t::ordinal_type &numCols,
      typename sparseMat_t::size_type &nnz,
      const int whichMatrix,
      typename sparseMat_t::ordinal_type &blockDim)
  {

    typedef typename sparseMat_t::StaticCrsGraphType::row_map_type::non_const_type ptr_type ;
    typedef typename sparseMat_t::StaticCrsGraphType::entries_type::non_const_type ind_type ;
    typedef typename sparseMat_t::values_type::non_const_type val_type ;
    typedef typename sparseMat_t::ordinal_type lno_t;
    typedef typename sparseMat_t::size_type size_type;
    typedef typename sparseMat_t::value_type scalar_t;

    using Kokkos::HostSpace;
    using Kokkos::MemoryUnmanaged;
    using Kokkos::View;

    if (whichMatrix == 0) {
      numRows = 8;
      numCols = 8;
      nnz = 24;
      blockDim = 1;

      const size_type ptrRaw[] = {0, 4, 8, 10, 12, 14, 16, 20, 24};
      const lno_t indRaw[] = {0, 1, 4, 5, 0, 1, 4, 5, 2, 3, 2, 3, 4, 5, 4, 5, 2, 3, 6, 7, 2, 3, 6, 7};
      const scalar_t valRaw[] = {.1, 1, 4, 5, -.1, -1, -4, -5, 2, 3, -2, -3, 4, 5, -4, -5, 2, 3, 6, 7, -2, -3, -6, -7};

      // Create the output Views.
      ptr = ptr_type("ptr", numRows + 1);
      ind = ind_type("ind", nnz);
      val = val_type("val", nnz);

      // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
      typename ptr_type::HostMirror::const_type  ptrIn( ptrRaw , numRows+1 );
      typename ind_type::HostMirror::const_type  indIn( indRaw , nnz );
      typename val_type::HostMirror::const_type  valIn( valRaw , nnz );

      Kokkos::deep_copy (ptr, ptrIn);
      Kokkos::deep_copy (ind, indIn);
      Kokkos::deep_copy (val, valIn);
    }
    else if (whichMatrix == 1) {
      numRows = 4;
      numCols = 4;
      nnz = 24;

      blockDim = 2;
      const lno_t numBlocks = 6;

      const size_type ptrRaw[] = {0, 2, 3, 4, 6};
      const lno_t indRaw[] = {0, 2, 1, 2, 1, 3};
      const scalar_t valRaw[] = {.1, 1, 4, 5, -.1, -1, -4, -5, 2, 3, -2, -3, 4, 5, -4, -5, 2, 3, 6, 7, -2, -3, -6, -7};

      // Create the output Views.
      ptr = ptr_type("ptr", numRows + 1);
      ind = ind_type("ind", numBlocks);
      val = val_type("val", nnz);

      // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
      typename ptr_type::HostMirror::const_type  ptrIn( ptrRaw , numRows+1 );
      typename ind_type::HostMirror::const_type  indIn( indRaw , numBlocks );
      typename val_type::HostMirror::const_type  valIn( valRaw , nnz );

      Kokkos::deep_copy (ptr, ptrIn);
      Kokkos::deep_copy (ind, indIn);
      Kokkos::deep_copy (val, valIn);
    }

    else { // whichMatrix != 0
      std::ostringstream os;
      os << "Invalid whichMatrix value " << whichMatrix
         << ".  Valid value(s) include " << 0 << ".";
      throw std::invalid_argument (os.str ());
    }
  }

  // Return the Kokkos::CrsMatrix corresponding to makeSparseMatrix().
  template<typename crsMat_t>
  crsMat_t  makeCrsMatrix_BlockStructure ()
  {
    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
    typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
    typedef typename crsMat_t::ordinal_type lno_t;
    typedef typename crsMat_t::size_type size_type;

    lno_view_t ptr;
    lno_nnz_view_t ind;
    scalar_view_t val;
    lno_t numRows;
    lno_t numCols;
    size_type nnz;
    lno_t blockDim;

    const int whichMatrix = 0;
    makeSparseMatrix<crsMat_t> (ptr, ind, val, numRows, numCols, nnz, whichMatrix, blockDim);
    return crsMat_t ("A", numRows, numCols, nnz, val, ptr, ind);
  }

  template<typename blkcrsMat_t>
  blkcrsMat_t  makeBlockCrsMatrix ()
  {
    typedef typename blkcrsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
    typedef typename blkcrsMat_t::values_type::non_const_type scalar_view_t;
    typedef typename blkcrsMat_t::ordinal_type lno_t;
    typedef typename blkcrsMat_t::size_type size_type;

    lno_view_t ptr;
    lno_nnz_view_t ind;
    scalar_view_t val;
    lno_t numRows;
    lno_t numCols;
    size_type nnz;
    lno_t blockDim;

    const int whichMatrix = 1;
    makeSparseMatrix<blkcrsMat_t> (ptr, ind, val, numRows, numCols, nnz, whichMatrix, blockDim);
    return blkcrsMat_t ("blkA", numRows, numCols, nnz, val, ptr, ind, blockDim);
  }


  template < class MatrixType, class ResultsType >
  struct TestFunctor {

  typedef typename MatrixType::value_type scalar_t;
  typedef typename MatrixType::ordinal_type lno_t;

  // Members
  MatrixType A;
  ResultsType d_results;

  // Constructor
  TestFunctor( MatrixType & A_, ResultsType & d_results_ ) :
    A(A_),
    d_results(d_results_)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int rid ) const
  {
    // Test 1: Check member functions behave as expected
    bool check0 = true;
    bool check1 = true;
    bool check2 = true;
    bool check3 = true;
    for ( lno_t i = 0; i < A.numRows(); ++i ) {

      // Test SparseBlockRowView
      {
        auto iblockrow = A.block_row(i);
        auto num_blocks_in_row = iblockrow.length;
        for ( auto blk = 0; blk < num_blocks_in_row; ++blk ){
          auto view_blk = iblockrow.block(blk);
          for ( auto lrow = 0; lrow < A.blockDim(); ++lrow ) {
            auto row_ptr = iblockrow.local_row_in_block(blk, lrow);
            for ( auto lcol = 0; lcol < A.blockDim(); ++lcol ) {
              auto entry = iblockrow.local_block_value(blk, lrow, lcol);
              //std::cout << "check0: " << ( entry == row_ptr[lcol] );
              //std::cout << "check1: " << ( entry == view_blk(lrow,lcol) );
              check0 = check0 && ( entry == row_ptr[lcol] );
              check1 = check1 && ( entry == view_blk(lrow,lcol) );
            } // end local col in row
          } // end local row in blk
        } // end blk
      }
      d_results(0) = check0;
      d_results(1) = check1;

      // Test SparseBlockRowViewConst
      {
        auto iblockrow = A.block_row_Const(i);
        auto num_blocks_in_row = iblockrow.length;
        for ( auto blk = 0; blk < num_blocks_in_row; ++blk ){
          auto view_blk = iblockrow.block(blk);
          for ( auto lrow = 0; lrow < A.blockDim(); ++lrow ) {
            auto row_ptr = iblockrow.local_row_in_block(blk, lrow);
            for ( auto lcol = 0; lcol < A.blockDim(); ++lcol ) {
              auto entry = iblockrow.local_block_value(blk, lrow, lcol);
              check2 = check2 && ( entry == row_ptr[lcol] );
              check3 = check3 && ( entry == view_blk(lrow,lcol) );
            } // end local col in row
          } // end local row in blk
        } // end blk
      }
      d_results(0) = check0;
      d_results(1) = check1;
      d_results(2) = check2;
      d_results(3) = check3;
    } // end for blk rows

    // Test sumIntoValues
    {
      check0 = true;
      check1 = true;
      check2 = true;
      const lno_t ncols = 1;
      const lno_t cols[] = {3};
      const lno_t browi = 3;
      const scalar_t vals[] = {10, 11, 20, 22}; // represents a single block: [10 11; 20 22]
      const scalar_t result[] = {16, 18, 14, 15};

      // This block will be summed into the existing block [6 7; -6 -7]
      // Expected result: [16 18; 14 15]
      A.sumIntoValues( browi, cols, ncols, vals );
      auto iblockrow = A.block_row_Const(browi);
      auto relBlk = iblockrow.findRelBlockOffset(cols[0]);
      auto view_blk = iblockrow.block(relBlk);
      for ( auto lrow = 0; lrow < A.blockDim(); ++lrow ) {
        auto row_ptr = iblockrow.local_row_in_block(relBlk, lrow);
        for ( auto lcol = 0; lcol < A.blockDim(); ++lcol ) {
          auto entry = iblockrow.local_block_value(relBlk, lrow, lcol);
          check0 = check0 && ( entry == row_ptr[lcol] );
          check1 = check1 && ( entry == view_blk(lrow,lcol) );
          check2 = check2 && ( entry == result[ lrow*A.blockDim() + lcol ] );
        } // end local col in row
      } // end local row in blk
      d_results(4) = check0;
      d_results(5) = check1;
      d_results(6) = check2;
    }

    // Test replaceValues
    {
      check0 = true;
      check1 = true;
      check2 = true;
      const lno_t ncols = 1;
      const lno_t cols[] = {3};
      const lno_t browi = 3;
      const scalar_t valsreplace[] = {-10, -11, -20, -22}; // represents a single block: [10 11; 20 22]

      // The existing block to be replaced was: [6 7; -6 -7]
      A.replaceValues( browi, cols, ncols, valsreplace );

      auto iblockrow = A.block_row_Const(browi);
      auto relBlk = iblockrow.findRelBlockOffset(cols[0]);
      auto view_blk = iblockrow.block(relBlk);
      for ( auto lrow = 0; lrow < A.blockDim(); ++lrow ) {
        auto row_ptr = iblockrow.local_row_in_block(relBlk, lrow);
        for ( auto lcol = 0; lcol < A.blockDim(); ++lcol ) {
          auto entry = iblockrow.local_block_value(relBlk, lrow, lcol);
          check0 = check0 && ( entry == row_ptr[lcol] );
          check1 = check1 && ( entry == view_blk(lrow,lcol) );
          check2 = check2 && ( entry == valsreplace[ lrow*A.blockDim() + lcol ] );
        } // end local col in row
      } // end local row in blk
      d_results(7) = check0;
      d_results(8) = check1;
      d_results(9) = check2;
    }

  }// end operator()(i)
  }; // end TestFunctor

} // namespace (anonymous)

// Create a CrsMatrix and BlockCrsMatrix and test member functions.
template <typename scalar_t, typename lno_t, typename size_type, typename device>
void
testBlockCrsMatrix ()
{
  using namespace Test;

  typedef KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crs_matrix_type;
  typedef KokkosSparse::Experimental::BlockCrsMatrix<scalar_t, lno_t, device, void, size_type> block_crs_matrix_type;

  crs_matrix_type crsA = makeCrsMatrix_BlockStructure<crs_matrix_type> ();
  block_crs_matrix_type A = makeBlockCrsMatrix<block_crs_matrix_type> ();

  const int num_entries = 10;
  typedef Kokkos::View< bool[num_entries], device > result_view_type;
  result_view_type d_results("d_results");
  auto h_results = Kokkos::create_mirror_view( d_results );

  Kokkos::parallel_for( "KokkosSparse::Test::BlockCrsMatrix", Kokkos::RangePolicy<typename device::execution_space>(0, 1), Test::TestFunctor< block_crs_matrix_type, result_view_type>( A, d_results ) );

  Kokkos::deep_copy( h_results, d_results );

  for ( decltype(h_results.extent(0)) i = 0; i < h_results.extent(0); ++i ) {
    EXPECT_EQ( h_results[i], true );
  }

}

#define EXECUTE_BLOCKCRS_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, sparse ## _ ## blkcrsmatrix ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  testBlockCrsMatrix<SCALAR, ORDINAL, OFFSET, DEVICE> (); \
}


#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_BLOCKCRS_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

#undef EXECUTE_BLOCKCRS_TEST

