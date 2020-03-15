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
#include "Tpetra_Details_CuSparseHandle_fwd.hpp"
#include "Tpetra_Details_CuSparseVector_fwd.hpp"
// I'm actually using the class' methods in this test, so I need to
// include its declaration header.
#include "Tpetra_Details_CuSparseMatrix.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosBlas1_scal.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

#if defined(KOKKOS_ENABLE_CUDA) && defined(HAVE_TPETRACORE_CUSPARSE)

  using crs_graph_type = Kokkos::StaticCrsGraph<
    int,  // local column index type
    Kokkos::LayoutLeft,
    Kokkos::Cuda,
    void, // MemoryTraits
    int   // row offset type
  >;

  using host_crs_graph_type = Kokkos::StaticCrsGraph<
    int,  // local column index type
    Kokkos::LayoutLeft,
    crs_graph_type::row_map_type::HostMirror::device_type,
    void, // MemoryTraits
    int   // row offset type
  >;

  host_crs_graph_type
  crsGraphHostMirrorView(const crs_graph_type& G)
  {
    auto ind_h = Kokkos::create_mirror_view(G.entries);
    Kokkos::deep_copy(ind_h, G.entries);

    // row_map is const, so we must make a deep copy.
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    crs_graph_type::row_map_type::non_const_type ptr_h
      (view_alloc("ptr_h", WithoutInitializing), G.row_map.extent(0));
    Kokkos::deep_copy(ptr_h, G.row_map);

    return host_crs_graph_type(ind_h, ptr_h);
  }

  template<class Scalar>
  using crs_matrix_type = KokkosSparse::CrsMatrix<
    Scalar,
    int,  // local column index type
    Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
    void, // MemoryTraits
    int   // row offset type
  >;

  template<class Scalar>
  using host_crs_matrix_type = KokkosSparse::CrsMatrix<
    Scalar,
    int,  // local column index type
    typename crs_matrix_type<Scalar>::values_type::HostMirror::device_type,
    void, // MemoryTraits
    int   // row offset type
  >;

  template<class Scalar>
  host_crs_matrix_type<Scalar>
  crsMatrixHostMirrorView(const crs_matrix_type<Scalar>& A)
  {
    auto G_h = crsGraphHostMirrorView(A.graph);
    auto val_h = Kokkos::create_mirror_view(A.values);
    Kokkos::deep_copy(val_h, A.values);
    return host_crs_matrix_type<Scalar>("A_h", A.numCols(), val_h, G_h);
  }

  template<class Scalar>
  using vector_type = Kokkos::View<
    Scalar*,
    Kokkos::LayoutLeft,
    Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>
  >;

  template<class Scalar>
  void
  referenceSparseMatVec(const Teuchos::ETransp op,
                        const Scalar alpha,
                        const crs_matrix_type<Scalar>& A,
                        const vector_type<const Scalar>& x,
                        const Scalar beta,
                        const vector_type<Scalar>& y)
  {
    auto x_h = Kokkos::create_mirror_view(x);
    Kokkos::deep_copy(x_h, x);
    auto y_h = Kokkos::create_mirror_view(y);
    Kokkos::deep_copy(y_h, y);
    auto A_h = crsMatrixHostMirrorView(A);

    if (alpha == Scalar{}) {
      if (beta == Scalar{}) {
        Kokkos::deep_copy(y_h, Scalar{});
      }
      else {
        KokkosBlas::scal(y_h, beta, y_h);
      }
    }
    else {
      const int numRows = A.numRows();
      for (int row = 0; row < numRows; ++row) {
        auto A_row = A_h.rowConst(row);
        Scalar sum {};
        for (int k = 0; k < A_row.length; ++k) {
          sum += A_row.value(k) * x_h(k);
        }
        if (beta == Scalar{}) {
          y_h[row] = sum;
        }
        else {
          y_h[row] = beta * y_h[row] + sum;
        }
      }
    }
    Kokkos::deep_copy(y, y_h);
  }

#endif

  void
  testCuSparseMatrixResult(bool& success, Teuchos::FancyOStream& out)
  {
#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRACORE_CUSPARSE)
    out << "Running this test requires enabling CUDA in Kokkos, "
      "and enabling the CUSPARSE TPL in Tpetra." << std::endl;
    TEUCHOS_ASSERT( false );
#else
    //using Tpetra::Details::CuSparseHandle;
    using Tpetra::Details::CuSparseMatrix;
    using Tpetra::Details::CuSparseMatrixVectorMultiplyAlgorithm;
    using Tpetra::Details::getCuSparseHandle;
    using Tpetra::Details::getCuSparseMatrix;
    using Tpetra::Details::getCuSparseVector;
    using Tpetra::Details::cuSparseMatrixVectorMultiply;
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    using std::endl;

    out << "Test CuSparseMatrix (actual sparse mat-vec results)" << endl;
    Teuchos::OSTab tab1(out);

    // Don't do cuSPARSE things unless a cuSPARSE handle is active.
    out << "Get cuSPARSE handle" << endl;
    Kokkos::Cuda execSpace;
    auto cuSparseHandle = getCuSparseHandle(execSpace);
    TEST_ASSERT( cuSparseHandle.get() != nullptr );
    if (! success) {
      return;
    }

    const int numRows = 3;
    const int numCols = 4;
    const int numEnt = 8;
    out << "Create nonsymmetric test matrix: numRows=" << numRows
        << ", numCols=" << numCols << ", numEnt=" << numEnt << endl;

    crs_graph_type::entries_type ind
      (view_alloc("ind", WithoutInitializing), numEnt);
    crs_matrix_type<double>::values_type val
      (view_alloc("val", WithoutInitializing), numEnt);
    crs_graph_type::row_map_type::non_const_type ptr
      (view_alloc("ptr", WithoutInitializing), numRows+1);

    auto val_h = Kokkos::create_mirror_view(val);
    auto ind_h = Kokkos::create_mirror_view(ind);
    auto ptr_h = Kokkos::create_mirror_view(ptr);
    {
      int pos = 0;
      int row = 0;
      ptr[row] = pos;
      // row 0

      val[pos] = 4.0;
      ind[pos] = 0;
      ++pos;

      val[pos] = 1.0;
      ind[pos] = 1;
      ++pos;

      ptr[row+1] = pos;
      ++row;
      // row 1

      val[pos] = -1.0;
      ind[pos] = 0;
      ++pos;

      val[pos] = 4.0;
      ind[pos] = 1;
      ++pos;

      val[pos] = 1.0;
      ind[pos] = 2;
      ++pos;

      ptr[row+1] = pos;
      ++row;
      // row 2

      val[pos] = -1.0;
      ind[pos] = 1;
      ++pos;

      val[pos] = 4.0;
      ind[pos] = 2;
      ++pos;

      val[pos] = 1.0;
      ind[pos] = 3;
      ++pos;

      ptr[row+1] = pos;
      ++row;

      // test
      TEUCHOS_ASSERT( pos == numEnt );
      TEUCHOS_ASSERT( row == numRows );

      Kokkos::deep_copy(val, val_h);
      Kokkos::deep_copy(ind, ind_h);
      Kokkos::deep_copy(ptr, ptr_h);
    }

    crs_graph_type G(ind, crs_graph_type::row_map_type(ptr));
    crs_matrix_type<double> A("A", numCols, val, G);

    vector_type<double> x(view_alloc("x", WithoutInitializing), numCols);
    vector_type<double> y(view_alloc("y", WithoutInitializing), numRows);
    auto y_h = Kokkos::create_mirror_view(y);

    const char* algStrings[] = {"DEFAULT", "LOAD_BALANCED"};
    for (const auto alg :
           {CuSparseMatrixVectorMultiplyAlgorithm::DEFAULT,
            CuSparseMatrixVectorMultiplyAlgorithm::LOAD_BALANCED}) {
      out << "1st pass: Test cuSparseMatrixVectorMultiply with alg="
          << algStrings[static_cast<int>(alg)] << endl;
      Kokkos::deep_copy(x, 1.0);
      auto A_cu = getCuSparseMatrix(numRows, numCols, numEnt,
                                    ptr.data(), ind.data(),
                                    val.data(), alg);
      auto x_cu = getCuSparseVector(x.data(), numCols);
      auto y_cu = getCuSparseVector(y.data(), numRows);

      const double alpha = 1.0;
      const double beta = 0.0;
      cuSparseMatrixVectorMultiply(*cuSparseHandle, Teuchos::NO_TRANS,
                                   alpha, *A_cu, *x_cu, beta, *y_cu);
      Kokkos::deep_copy(y_h, y);
      TEST_EQUALITY( y_h[0], 5.0 );
      TEST_EQUALITY( y_h[1], 4.0 );
      TEST_EQUALITY( y_h[2], 4.0 );
    }

    out << "Test whether we can change the matrix's values and "
      "still get the right answer" << endl;

    Kokkos::deep_copy(val_h, val);
    // Change diagonal entry of row 1.
    val_h[3] = 6.0;
    Kokkos::deep_copy(val, val_h);

    for (const auto alg :
           {CuSparseMatrixVectorMultiplyAlgorithm::DEFAULT,
            CuSparseMatrixVectorMultiplyAlgorithm::LOAD_BALANCED}) {
      out << "2nd pass: Test cuSparseMatrixVectorMultiply with alg="
          << algStrings[static_cast<int>(alg)] << endl;
      Kokkos::deep_copy(x, 1.0);
      auto A_cu = getCuSparseMatrix(numRows, numCols, numEnt,
                                    ptr.data(), ind.data(),
                                    val.data(), alg);
      auto x_cu = getCuSparseVector(x.data(), numCols);
      auto y_cu = getCuSparseVector(y.data(), numRows);

      const double alpha = 1.0;
      const double beta = 0.0;
      cuSparseMatrixVectorMultiply(*cuSparseHandle, Teuchos::NO_TRANS,
                                   alpha, *A_cu, *x_cu, beta, *y_cu);
      Kokkos::deep_copy(y_h, y);
      TEST_EQUALITY( y_h[0], 5.0 );
      TEST_EQUALITY( y_h[1], 6.0 );
      TEST_EQUALITY( y_h[2], 4.0 );
    }

    out << "Done with test!" << endl;
#endif
  }

  void
  testCuSparseMatrix(bool& success, Teuchos::FancyOStream& out)
  {
#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRACORE_CUSPARSE)
    out << "Running this test requires enabling CUDA in Kokkos, "
      "and enabling the CUSPARSE TPL in Tpetra." << std::endl;
    TEUCHOS_ASSERT( false );
#else
    using Tpetra::Details::CuSparseHandle;
    using Tpetra::Details::CuSparseMatrix;
    using Tpetra::Details::CuSparseMatrixVectorMultiplyAlgorithm;
    using Tpetra::Details::getCuSparseMatrix;
    using Tpetra::Details::getCuSparseVector;
    using Tpetra::Details::cuSparseMatrixVectorMultiply;
    // using Kokkos::view_alloc;
    // using Kokkos::WithoutInitializing;
    using std::endl;
    using LO = Tpetra::Details::DefaultTypes::local_ordinal_type;

    out << "Test CuSparseMatrix" << endl;
    Teuchos::OSTab tab1(out);

    // Don't do cuSPARSE things unless a cuSPARSE handle is active.
    Kokkos::Cuda execSpace1;
    std::shared_ptr<CuSparseHandle> h1 =
      Tpetra::Details::getCuSparseHandle(execSpace1);
    TEST_ASSERT( h1.get() != nullptr );
    if (! success) {
      return;
    }

    using memory_space = Kokkos::CudaUVMSpace;

    for (const auto alg :
           {CuSparseMatrixVectorMultiplyAlgorithm::DEFAULT,
            CuSparseMatrixVectorMultiplyAlgorithm::LOAD_BALANCED}) {
      out << "Test alg="
          << (alg == CuSparseMatrixVectorMultiplyAlgorithm::DEFAULT
              ? "DEFAULT" : "LOAD_BALANCED") << endl;
      Teuchos::OSTab tab1_1(out);

      for (int64_t numRows : {0, 1, 11}) {
        out << "numRows: " << numRows << endl;
        Teuchos::OSTab tab2(out);

        Kokkos::View<float*, memory_space> y_f_k("y_f", numRows);
        auto y_f = getCuSparseVector(y_f_k.data(), numRows);

        Kokkos::View<double*, memory_space> y_d_k("y_d", numRows);
        auto y_d = getCuSparseVector(y_d_k.data(), numRows);

        Kokkos::View<LO*, memory_space> ptr("ptr", numRows+1);
        for(int64_t numEntries : {0, 32}) {
          out << "numEntries: " << numEntries << endl;
          Teuchos::OSTab tab3(out);

          Kokkos::View<LO*, memory_space> ind("ind", numEntries);
          Kokkos::View<float*, memory_space> val_f
            ("val_f", numEntries);
          Kokkos::View<double*, memory_space> val_d
            ("val_d", numEntries);

          for (int64_t numCols : {0, 1, 5, 13}) {
            const int64_t numEnt = (numRows == 0 || numCols == 0) ?
              int64_t(0) : numEntries;
            out << "numCols: " << numCols << ", numEnt: " << numEnt
                << endl;
            Teuchos::OSTab tab4(out);

            Kokkos::View<float*, memory_space> x_f_k("x_f", numCols);
            auto x_f = getCuSparseVector(x_f_k.data(), numCols);

            Kokkos::View<double*, memory_space> x_d_k("x_d", numCols);
            auto x_d = getCuSparseVector(x_d_k.data(), numCols);

            out << "Call getCuSparseMatrix (float)" << endl;
            auto mat_f = getCuSparseMatrix(numRows, numCols, numEnt,
                                           ptr.data(), ind.data(),
                                           val_f.data(), alg);
            out << "Test result of getCuSparseMatrix (float)" << endl;
            TEST_ASSERT( mat_f.get() != nullptr );
            if (mat_f.get() != nullptr) {
              cusparseMatDescr_t descr_f = mat_f->getDescr();
              out << "mat_f->getDescr() returned" << endl;
              TEST_ASSERT( cusparseGetMatType(descr_f) ==
                           CUSPARSE_MATRIX_TYPE_GENERAL );
              TEST_ASSERT( cusparseGetMatDiagType(descr_f) ==
                           CUSPARSE_DIAG_TYPE_NON_UNIT );
              TEST_ASSERT( cusparseGetMatIndexBase(descr_f) ==
                           CUSPARSE_INDEX_BASE_ZERO );

              out << "Call cuSparseMatrixVectorMultiply (float)" << endl;
              const float alpha_f (1.2);
              const float beta_f (2.3);
              cuSparseMatrixVectorMultiply(*h1, Teuchos::NO_TRANS,
                                           alpha_f, *mat_f, *x_f,
                                           beta_f, *y_f);
            }

            out << "Call getCuSparseMatrix (double)" << endl;
            auto mat_d = getCuSparseMatrix(numRows, numCols, numEnt,
                                           ptr.data(), ind.data(),
                                           val_d.data(), alg);
            out << "Test result of getCuSparseMatrix (double)" << endl;
            TEST_ASSERT( mat_d.get() != nullptr );
            if (mat_d.get() != nullptr) {
              cusparseMatDescr_t descr_d = mat_d->getDescr();
              out << "mat_d->getDescr() returned" << endl;
              TEST_ASSERT( cusparseGetMatType(descr_d) ==
                           CUSPARSE_MATRIX_TYPE_GENERAL );
              TEST_ASSERT( cusparseGetMatDiagType(descr_d) ==
                           CUSPARSE_DIAG_TYPE_NON_UNIT );
              TEST_ASSERT( cusparseGetMatIndexBase(descr_d) ==
                           CUSPARSE_INDEX_BASE_ZERO );

              out << "Call cuSparseMatrixVectorMultiply (double)" << endl;
              const float alpha_d (1.2);
              const float beta_d (2.3);
              cuSparseMatrixVectorMultiply(*h1, Teuchos::NO_TRANS,
                                           alpha_d, *mat_d, *x_d,
                                           beta_d, *y_d);
            }
          }
        }
      }
    }
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

  TEUCHOS_UNIT_TEST( Utils, CuSparseMatrix )
  {
    runTest(success, out, testCuSparseMatrix);
  }

  TEUCHOS_UNIT_TEST( Utils, CuSparseMatrixResult )
  {
    runTest(success, out, testCuSparseMatrixResult);
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
