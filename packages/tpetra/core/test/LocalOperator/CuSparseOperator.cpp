/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_CuSparseHandle_fwd.hpp"
#include "Tpetra_Details_LocalCuSparseCrsMatrixOperator.hpp"
#include "Tpetra_Details_LocalCrsMatrixOperatorWithSetup.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse.hpp"
#include "Kokkos_Core.hpp"
#include <type_traits>

namespace { // (anonymous)

#if defined(KOKKOS_ENABLE_CUDA) && defined(HAVE_TPETRACORE_CUSPARSE)
  using device_type =
    Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>;

  template<class Scalar>
  using crs_matrix_type =
    typename Tpetra::Details::LocalCrsMatrixOperatorWithSetup<
      Scalar, Scalar, device_type>::local_matrix_type;

  using local_ordinal_type =
    crs_matrix_type<double>::staticcrsgraph_type::data_type;

  template<class Scalar>
  using host_crs_matrix_type = KokkosSparse::CrsMatrix<
    Scalar,
    local_ordinal_type,
    typename crs_matrix_type<Scalar>::values_type::HostMirror::device_type,
    void,   // MemoryTraits
    typename crs_matrix_type<Scalar>::size_type
  >;

  using crs_graph_type =
    typename crs_matrix_type<double>::staticcrsgraph_type;

  using host_crs_graph_type =
    typename host_crs_matrix_type<double>::staticcrsgraph_type;

  template<class Scalar>
  using multivector_type = Kokkos::View<
    Scalar**,
    Kokkos::LayoutLeft,
    Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>
  >;

#endif

  void
  testCuSparseOperatorResult(bool& success, Teuchos::FancyOStream& out)
  {
#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRACORE_CUSPARSE)
    out << "Running this test requires enabling CUDA in Kokkos, "
      "and enabling the CUSPARSE TPL in Tpetra." << std::endl;
    TEUCHOS_ASSERT( false );
#else
    using Tpetra::Details::getCuSparseHandle;
    using Tpetra::Details::LocalCuSparseCrsMatrixOperator;
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    using std::endl;
    using row_offset_type = crs_graph_type::size_type;

    out << "Test LocalCuSparseCrsMatrixOperator with row_offset_type="
        << Teuchos::TypeNameTraits<row_offset_type>::name() << endl;
    Teuchos::OSTab tab1(out);

    // Don't do cuSPARSE things unless a cuSPARSE handle is active.
    out << "Get cuSPARSE handle" << endl;
    Kokkos::Cuda execSpace;
    auto cuSparseHandle = getCuSparseHandle(execSpace);
    TEST_ASSERT( cuSparseHandle.get() != nullptr );
    if (! success) {
      return;
    }

    const local_ordinal_type numRows = 3;
    const local_ordinal_type numCols = 4;
    const local_ordinal_type numEnt = 8;
    out << "Create nonsymmetric test matrix: numRows=" << numRows
        << ", numCols=" << numCols << ", numEnt=" << numEnt << endl;

    using inds_type = crs_graph_type::entries_type;
    using vals_type = crs_matrix_type<double>::values_type;
    using offsets_type = crs_graph_type::row_map_type::non_const_type;

    inds_type ind(view_alloc("ind", WithoutInitializing), numEnt);
    vals_type val(view_alloc("val", WithoutInitializing), numEnt);
    offsets_type ptr(view_alloc("ptr", WithoutInitializing), numRows+1);

    auto val_h = Kokkos::create_mirror_view(val);
    auto ind_h = Kokkos::create_mirror_view(ind);
    auto ptr_h = Kokkos::create_mirror_view(ptr);
    {
      row_offset_type pos = 0;
      local_ordinal_type row = 0;
      ptr_h[row] = pos;
      // row 0

      val_h[pos] = 4.0;
      ind_h[pos] = 0;
      ++pos;

      val_h[pos] = 1.0;
      ind_h[pos] = 1;
      ++pos;

      ptr_h[row+1] = pos;
      ++row;
      // row 1

      val_h[pos] = -1.0;
      ind_h[pos] = 0;
      ++pos;

      val_h[pos] = 4.0;
      ind_h[pos] = 1;
      ++pos;

      val_h[pos] = 1.0;
      ind_h[pos] = 2;
      ++pos;

      ptr_h[row+1] = pos;
      ++row;
      // row 2

      val_h[pos] = -1.0;
      ind_h[pos] = 1;
      ++pos;

      val_h[pos] = 4.0;
      ind_h[pos] = 2;
      ++pos;

      val_h[pos] = 1.0;
      ind_h[pos] = 3;
      ++pos;

      ptr_h[row+1] = pos;
      ++row;

      // test
      TEST_ASSERT( pos == numEnt );
      TEST_ASSERT( row == numRows );

      {
        std::vector<double> x(numCols, 1.0);
        std::vector<double> y(numRows);

        for (int row = 0; row < numRows; ++row) {
          double sum = 0;
          for (int k = static_cast<int>(ptr_h[row]);
               k < static_cast<int>(ptr_h[row+1]); ++k) {
            sum += val_h[k] * x[ind_h[k]];
          }
          y[row] = sum;
        }

        TEST_EQUALITY( y[0], 5.0 );
        TEST_EQUALITY( y[1], 4.0 );
        TEST_EQUALITY( y[2], 4.0 );

        if (! success) {
          out << "Test not successful thus far; exiting early." << endl;
          return;
        }
      }

      Kokkos::deep_copy(val, val_h);
      Kokkos::deep_copy(ind, ind_h);
      Kokkos::deep_copy(ptr, ptr_h);
    }

    crs_graph_type G(ind, crs_graph_type::row_map_type(ptr));

    static_assert(
      std::is_same<crs_graph_type,
        typename crs_matrix_type<double>::staticcrsgraph_type
      >::value, "StaticCrsGraph type mismatch.");

    using mat_type = crs_matrix_type<double>;
    std::shared_ptr<mat_type> A(new mat_type(std::string("A"), numCols, val, G));

    TEST_EQUALITY( A->numRows(), numRows );
    TEST_EQUALITY( A->numCols(), numCols );

    out << "Before creating the LocalCuSparseCrsMatrixOperator, "
      "test that KokkosSparse::spmv gives the right answer." << endl;
    for (const local_ordinal_type numVecs : {1, 7}) {
      Teuchos::OSTab tab2(out);
      out << "numVecs=" << numVecs << endl;

      multivector_type<double> X
        (view_alloc("X", WithoutInitializing), numCols, numVecs);
      multivector_type<double> Y
        (view_alloc("Y", WithoutInitializing), numRows, numVecs);
      auto Y_h = Kokkos::create_mirror_view(Y);

      {
        Teuchos::OSTab tab3(out);
        out << "Compare to KokkosSparse::spmv" << endl;
        Kokkos::deep_copy(X, 1.0);
        const double alpha = 1.0;
        const double beta = 0.0;
        KokkosSparse::spmv(KokkosSparse::NoTranspose,
                           alpha, *A, X, beta, Y);
        Kokkos::deep_copy(Y_h, Y);
        for (int vec = 0; vec < numVecs; ++vec) {
          TEST_EQUALITY( Y_h(0, vec), 5.0 );
          TEST_EQUALITY( Y_h(1, vec), 4.0 );
          TEST_EQUALITY( Y_h(2, vec), 4.0 );
        }
      }

      if (! success) {
        out << "Test not successful thus far; exiting early." << endl;
        return;
      }
    }

    LocalCuSparseCrsMatrixOperator<double> A_op(execSpace, A);
    TEST_ASSERT( ! A_op.isFillComplete() );

    TEST_NOTHROW( A_op.fillComplete() );
    TEST_ASSERT( A_op.isFillComplete() );

    out << "Test original matrix, after operator creation" << endl;

    for (const local_ordinal_type numVecs : {1, 7}) {
      Teuchos::OSTab tab2(out);
      out << "numVecs=" << numVecs << endl;

      multivector_type<double> X
        (view_alloc("X", WithoutInitializing), numCols, numVecs);
      multivector_type<double> Y
        (view_alloc("Y", WithoutInitializing), numRows, numVecs);
      auto Y_h = Kokkos::create_mirror_view(Y);

      {
        Teuchos::OSTab tab3(out);
        out << "Compare to KokkosSparse::spmv" << endl;
        Kokkos::deep_copy(X, 1.0);
        const double alpha = 1.0;
        const double beta = 0.0;
        KokkosSparse::spmv(KokkosSparse::NoTranspose,
                           alpha, *A, X, beta, Y);
        Kokkos::deep_copy(Y_h, Y);
        for (int vec = 0; vec < numVecs; ++vec) {
          TEST_EQUALITY( Y_h(0, vec), 5.0 );
          TEST_EQUALITY( Y_h(1, vec), 4.0 );
          TEST_EQUALITY( Y_h(2, vec), 4.0 );
        }
      }

      if (! success) {
        out << "Test not successful thus far; exiting early." << endl;
        return;
      }
      Kokkos::deep_copy(Y, -1.0);

      {
        Teuchos::OSTab tab3(out);
        out << "Test LocalCuSparseCrsMatrixOperator::apply" << endl;
        Kokkos::deep_copy(X, 1.0);

        {
          const double alpha = 1.0;
          const double beta = 0.0;
          out << "alpha=" << alpha << ", beta=" << beta << endl;
          Teuchos::OSTab tab4(out);

          A_op.apply(X, Y, Teuchos::NO_TRANS, alpha, beta);
          Kokkos::deep_copy(Y_h, Y);
          for (int vec = 0; vec < numVecs; ++vec) {
            TEST_EQUALITY( Y_h(0, vec), 5.0 );
            TEST_EQUALITY( Y_h(1, vec), 4.0 );
            TEST_EQUALITY( Y_h(2, vec), 4.0 );
          }
        }
        {
          const double alpha = 0.0;
          const double beta = 1.0;
          out << "alpha=" << alpha << ", beta=" << beta << endl;
          Teuchos::OSTab tab4(out);

          Kokkos::deep_copy(Y, 7.0);
          A_op.apply(X, Y, Teuchos::NO_TRANS, alpha, beta);
          Kokkos::deep_copy(Y_h, Y);
          for (int vec = 0; vec < numVecs; ++vec) {
            TEST_EQUALITY( Y_h(0, vec), 7.0 );
            TEST_EQUALITY( Y_h(1, vec), 7.0 );
            TEST_EQUALITY( Y_h(2, vec), 7.0 );
          }
        }
        {
          const double alpha = 2.0;
          const double beta = -1.0;
          out << "alpha=" << alpha << ", beta=" << beta << endl;
          Teuchos::OSTab tab4(out);

          Kokkos::deep_copy(Y, 7.0);
          A_op.apply(X, Y, Teuchos::NO_TRANS, alpha, beta);
          Kokkos::deep_copy(Y_h, Y);
          for (int vec = 0; vec < numVecs; ++vec) {
            TEST_EQUALITY( Y_h(0, vec), 3.0 );
            TEST_EQUALITY( Y_h(1, vec), 1.0 );
            TEST_EQUALITY( Y_h(2, vec), 1.0 );
          }
        }
      }

      if (! success) {
        out << "Test not successful thus far; exiting early." << endl;
        return;
      }
    }

    out << "Call A_op.resumeFill() and modify the matrix's values, "
      "without recreating the operator." << endl;

    TEST_NOTHROW( A_op.resumeFill() );
    TEST_ASSERT( ! A_op.isFillComplete() );

    Kokkos::deep_copy(val_h, val);
    // Change diagonal entry of row 1.
    val_h[3] = 6.0;
    Kokkos::deep_copy(val, val_h);

    out << "Call A_op.fillComplete()" << endl;

    TEST_NOTHROW( A_op.fillComplete() );
    TEST_ASSERT( A_op.isFillComplete() );

    for (const local_ordinal_type numVecs : {1, 7}) {
      Teuchos::OSTab tab2(out);
      out << "numVecs=" << numVecs << endl;

      multivector_type<double> X
        (view_alloc("X", WithoutInitializing), numCols, numVecs);
      multivector_type<double> Y
        (view_alloc("Y", WithoutInitializing), numRows, numVecs);
      auto Y_h = Kokkos::create_mirror_view(Y);

      {
        Teuchos::OSTab tab3(out);
        out << "Compare to KokkosSparse::spmv" << endl;
        Kokkos::deep_copy(X, 1.0);
        const double alpha = 1.0;
        const double beta = 0.0;
        KokkosSparse::spmv(KokkosSparse::NoTranspose,
                           alpha, *A, X, beta, Y);
        Kokkos::deep_copy(Y_h, Y);
        for (int vec = 0; vec < numVecs; ++vec) {
          TEST_EQUALITY( Y_h(0, vec), 5.0 );
          TEST_EQUALITY( Y_h(1, vec), 6.0 );
          TEST_EQUALITY( Y_h(2, vec), 4.0 );
        }
      }

      Kokkos::deep_copy(Y, -1.0);
      {
        Teuchos::OSTab tab3(out);
        out << "Test LocalCuSparseCrsMatrixOperator::apply" << endl;
        Kokkos::deep_copy(X, 1.0);

        {
          const double alpha = 1.0;
          const double beta = 0.0;
          out << "alpha=" << alpha << ", beta=" << beta << endl;
          Teuchos::OSTab tab4(out);

          A_op.apply(X, Y, Teuchos::NO_TRANS, alpha, beta);
          Kokkos::deep_copy(Y_h, Y);
          for (int vec = 0; vec < numVecs; ++vec) {
            TEST_EQUALITY( Y_h(0, vec), 5.0 );
            TEST_EQUALITY( Y_h(1, vec), 6.0 );
            TEST_EQUALITY( Y_h(2, vec), 4.0 );
          }
        }
        {
          const double alpha = 0.0;
          const double beta = 1.0;
          out << "alpha=" << alpha << ", beta=" << beta << endl;
          Teuchos::OSTab tab4(out);

          Kokkos::deep_copy(Y, 7.0);
          A_op.apply(X, Y, Teuchos::NO_TRANS, alpha, beta);
          Kokkos::deep_copy(Y_h, Y);
          for (int vec = 0; vec < numVecs; ++vec) {
            TEST_EQUALITY( Y_h(0, vec), 7.0 );
            TEST_EQUALITY( Y_h(1, vec), 7.0 );
            TEST_EQUALITY( Y_h(2, vec), 7.0 );
          }
        }
        {
          const double alpha = 2.0;
          const double beta = -1.0;
          out << "alpha=" << alpha << ", beta=" << beta << endl;
          Teuchos::OSTab tab4(out);

          Kokkos::deep_copy(Y, 7.0);
          A_op.apply(X, Y, Teuchos::NO_TRANS, alpha, beta);
          Kokkos::deep_copy(Y_h, Y);
          for (int vec = 0; vec < numVecs; ++vec) {
            TEST_EQUALITY( Y_h(0, vec), 3.0 );
            TEST_EQUALITY( Y_h(1, vec), 5.0 );
            TEST_EQUALITY( Y_h(2, vec), 1.0 );
          }
        }
      }
    }

    out << "Done with test!" << endl;
#endif
  }

  void
  runTest(bool& success,
          Teuchos::FancyOStream& out,
          std::function<void(bool&, Teuchos::FancyOStream&)> testFunction)
  {
    // Replace 'out' with cerr in verbose mode.  This lets us diagnose
    // crashes, since the Teuchos unit test framework normally
    // captures output until the test finishes.
    using Teuchos::FancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    RCP<FancyOStream> myOutPtr;
    const bool verbose = Tpetra::Details::Behavior::verbose();
    if (verbose) {
      myOutPtr = Teuchos::getFancyOStream(rcpFromRef(std::cerr));
    }
    else {
      myOutPtr = rcpFromRef(out);
    }
    FancyOStream& myOut = *myOutPtr;
    testFunction(success, myOut);
  }

  TEUCHOS_UNIT_TEST( Utils, LocalCuSparseCrsMatrixOperator_apply_resume_apply )
  {
    runTest(success, out, &testCuSparseOperatorResult);
  }

} // namespace (anonymous)

int
main(int argc, char* argv[])
{
  int errCode = 0;
  {
    Kokkos::ScopeGuard kokkosScope(argc, argv);
    errCode = Teuchos::UnitTestRepository::
      runUnitTestsFromMain(argc, argv);
  }
  return errCode;
}
