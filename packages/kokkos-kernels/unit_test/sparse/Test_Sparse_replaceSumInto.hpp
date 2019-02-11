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

//#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <gtest/gtest.h>
#include "KokkosSparse_CrsMatrix.hpp"

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

// mfh 21 Jun 2016: CUDA 7.5 with GCC 4.8.4 gives me funny build
// errors if I put this functor in an anonymous namespace.  If I name
// the namespace, it builds just fine.
namespace Test {
  template<class CrsMatrixType>
  class ModifyEvenNumberedRows {
  public:
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef typename CrsMatrixType::value_type value_type;

    ModifyEvenNumberedRows (const CrsMatrixType& A,
                            const bool replace,
                            const bool sorted,
                            const bool atomic) :
      A_ (A), replace_ (replace), sorted_ (sorted)
    {}

    KOKKOS_FUNCTION void
    operator () (const ordinal_type& lclRow) const
    {
      if (lclRow % static_cast<ordinal_type> (2) == 0) {
        const ordinal_type ncol = 1;
        ordinal_type cols[1];
        value_type vals[1];

        const value_type ONE = Kokkos::Details::ArithTraits<value_type>::one ();
        const value_type THREE = ONE + ONE + ONE;

        cols[0] = lclRow;
        vals[0] = replace_ ? THREE : ONE;

        if (replace_) {
          A_.replaceValues (lclRow, cols, ncol, vals, sorted_, atomic_);
        }
        else { // sumInto
          A_.sumIntoValues (lclRow, cols, ncol, vals, sorted_, atomic_);
        }
      }
    }

  private:
    CrsMatrixType A_;
    bool replace_;
    bool sorted_;
    bool atomic_;
  };
} // namespace KokkosSparseTest

namespace { // (anonymous)
  using std::endl;

  template<class CrsMatrixType>
  void
  modifyEvenNumberedRows (const CrsMatrixType& A,
                          const bool replace,
                          const bool sorted,
                          const bool atomic)
  {
    typedef typename CrsMatrixType::device_type::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, typename CrsMatrixType::ordinal_type> policy_type;

    ::Test::ModifyEvenNumberedRows<CrsMatrixType> functor (A, replace, sorted, atomic);
    Kokkos::parallel_for ( "KokkosSparse::Test::ReplaceSumInto", policy_type (0, A.numRows ()), functor);
  }

  template<class CrsMatrixType>
  bool
  checkWhetherEvenNumberedRowsWereModified (const CrsMatrixType& A,
                                            const bool replace,
                                            const bool /* sorted */,
                                            const bool /* atomic */ )
  {
    typedef typename CrsMatrixType::value_type SC;
    typedef typename CrsMatrixType::ordinal_type LO;

    const SC ONE = Kokkos::Details::ArithTraits<SC>::one ();
    const SC TWO = ONE + ONE;
    const SC THREE = ONE + ONE + ONE;

    typename CrsMatrixType::values_type val = A.values;
    typename CrsMatrixType::values_type::HostMirror val_h = Kokkos::create_mirror_view (val);
    Kokkos::deep_copy (val_h, val);
    Kokkos::fence();
    const LO numRows = A.numRows ();
    bool success = true;
    for (LO lclRow = 0; lclRow < numRows; ++lclRow) {
      if (lclRow % 2 == 0) {
        if (replace && val_h(lclRow) != THREE) {
          success = false;
          break;
        }
        else if (! replace && val_h(lclRow) != TWO) {
          success = false;
          break;
        }
      }
      else if (val_h(lclRow) != ONE) {
        success = false;
        break;
      }
    }
    return success;
  }

  template<class CrsMatrixType>
  void
  testOneCase (bool& success,
               //Teuchos::FancyOStream& out,
          	   std::ostream &out,
               const CrsMatrixType& A,
               const bool replace,
               const bool sorted,
               const bool atomic)
  {
    using Kokkos::Details::ArithTraits;
    typedef typename CrsMatrixType::value_type value_type;

    //Teuchos::OSTab tab0 (out);
    out << "replace: " << (replace ? "true" : "false")
        << ", sorted: " << (sorted ? "true" : "false")
        << ", atomic: " << (atomic ? "true" : "false")
        << endl;

    modifyEvenNumberedRows (A, replace, sorted, atomic);
    const bool lclSuccess =
      checkWhetherEvenNumberedRowsWereModified (A, replace, sorted, atomic);
    EXPECT_TRUE( lclSuccess ); // this modifies 'success' and prints to 'out'
    // Restore original values.
    Kokkos::deep_copy (A.values, ArithTraits<value_type>::one ());
  }

  // Test findRelOffset with various array data types and for various cases.
  //
  // This takes the same arguments as if it were declared via the
  // TEUCHOS_UNIT_TEST macro.
  template <typename scalar_t, typename lno_t, typename size_type, typename device>
  void generalTest (bool& success, std::ostream &out)
		  	  	  	  	  	  	  //Teuchos::FancyOStream& out)
  {

    typedef KokkosSparse::CrsMatrix<scalar_t, lno_t, device,void, size_type> matrix_type;

    //Teuchos::OSTab tab0 (out);
    out << "Test KokkosSparse::CrsMatrix::{replace,sumInto}Values*" << endl;
    //Teuchos::OSTab tab1 (out);

    out << "Create a diagonal matrix as a test problem" << endl;

    const lno_t numRows = 10;
    typename matrix_type::size_type numEnt = 0; // to be updated below
    typename matrix_type::row_map_type::non_const_type ptr ("ptr", numRows+1);
    {
      typename matrix_type::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
      ptr_h[0] = 0;
      for (lno_t lclRow = 0; lclRow < numRows; ++lclRow) {
        ptr_h[lclRow+1] = ptr_h[lclRow] + 1; // 1 entry in each row
      }
      numEnt = ptr_h[numRows];
      Kokkos::deep_copy (ptr, ptr_h);
    }

    typename matrix_type::index_type::non_const_type ind ("ind", numEnt);
    {
      typename matrix_type::index_type::HostMirror ind_h = Kokkos::create_mirror_view (ind);
      for (lno_t lclRow = 0; lclRow < numRows; ++lclRow) {
        ind_h[lclRow] = lclRow; // diagonal matrix
      }
      Kokkos::deep_copy (ind, ind_h);
    }

    typename matrix_type::values_type val ("val", numEnt);
    {
      typename matrix_type::values_type::HostMirror val_h = Kokkos::create_mirror_view (val);
      for (lno_t lclRow = 0; lclRow < numRows; ++lclRow) {
        val_h[lclRow] = 1.0; // diagonal matrix
      }
      Kokkos::deep_copy (val, val_h);
    }

    const lno_t numCols = numRows; // square diagonal matrix
    matrix_type A ("A", numRows, numCols, numEnt, val, ptr, ind);

    for (int replaceInt = 0; replaceInt < 2; ++replaceInt) {
      const bool replace = replaceInt != 0;
      for (int sortedInt = 0; sortedInt < 2; ++sortedInt) {
        const bool sorted = sortedInt != 0;
        for (int atomicInt = 0; atomicInt < 2; ++atomicInt) {
          const bool atomic = atomicInt != 0;
          testOneCase (success, out, A, replace, sorted, atomic);
        }
      }
    }
  }

} // namespace (anonymous)



template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_replaceSumInto()
{
  using std::endl;
  using std::endl;
  class NullBuffer : public std::streambuf
  {
  public:
    int overflow(int c) { return c; }
  };
  NullBuffer null_buffer;
  //std::ostream &out = std::cout;
  std::ostream out(&null_buffer);

  bool success = true;
  out << "Run test" << endl;
  generalTest <scalar_t, lno_t, size_type, device> (success, out);
  EXPECT_TRUE( success);
}


#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## replaceSumInto ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_replaceSumInto<SCALAR,ORDINAL,OFFSET,DEVICE>(); \
}


#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif





