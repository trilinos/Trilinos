#include <gtest/gtest.h>
#include "Kokkos_Core.hpp"
#include "KokkosSparse_spmv.hpp"
#include <limits>
#include <type_traits>

namespace { // (anonymous)

// See GitHub Issue #101 for details:
//
// https://github.com/kokkos/kokkos-kernels/issues/101
//
// Create a 1x2 matrix A = [1, eps_f/2], where eps_f is machine
// epsilon for float, and multiply it by x = [1, 1]^T.  If y = A*x and
// x and y are both in double precision, then KokkosSparse::spmv
// should compute the result in double precision ("double"), not
// single precision ("float").  If it computes in double precision,
// then the result should not equal 1 (or 0); if it computes in single
// precison, then the result should equal 1.
template<class DeviceType>
void
test_github_issue_101 ()
{
  typedef KokkosSparse::CrsMatrix<float, int, DeviceType> float_matrix_type;
  typedef KokkosSparse::CrsMatrix<double, int, DeviceType> double_matrix_type;
  static_assert (std::is_same<typename float_matrix_type::StaticCrsGraphType,
    typename double_matrix_type::StaticCrsGraphType>::value,
    "Two KokkosSparse::CrsMatrix types that differ only in the type of "
    "matrix values, appear to have two different StaticCrsGraphType "
    "typedefs.  This should never happen.");
  typedef typename float_matrix_type::StaticCrsGraphType graph_type;

  constexpr int numRows = 1;
  constexpr int numCols = 2;
  constexpr double alpha_d = 1.0;
  constexpr double beta_d = 0.0;
  const float EPS_f = std::numeric_limits<float>::epsilon ();

  graph_type G;
  {
    typename graph_type::entries_type colInds ("colInds", numCols);
    auto colInds_h = Kokkos::create_mirror_view (colInds);
    colInds_h[0] = 0;
    colInds_h[1] = 1;
    Kokkos::deep_copy (colInds, colInds_h);

    typedef typename graph_type::row_map_type::non_const_type row_offsets_type;
    row_offsets_type rowOffsets ("rowOffsets", numRows+1);
    auto rowOffsets_h = Kokkos::create_mirror_view (rowOffsets);
    rowOffsets_h[0] = 0; // Entries start at offset 0
    rowOffsets_h[1] = 2; // 2 entries total in the "sparse" matrix
    Kokkos::deep_copy (rowOffsets, rowOffsets_h);

    G = graph_type (colInds, rowOffsets);
  }

  Kokkos::View<double*, DeviceType> x ("x", numCols);
  Kokkos::deep_copy (x, static_cast<double> (1.0));
  Kokkos::View<double*, DeviceType> y ("y", numRows);
  auto y_h = Kokkos::create_mirror_view (y); // we'll want this later

  // Pick some number large enough to exercise all unrolling cases.
  // Sparse mat-vec does or at least used to unroll for 1, 2, ..., 17
  // vectors.  Include a little extra in case the implementers decide
  // to strip-mine that.
  constexpr int numVecs = 22;
  Kokkos::View<double**, Kokkos::LayoutLeft, DeviceType> X ("X", numCols, numVecs);
  Kokkos::deep_copy (X, static_cast<double> (1.0));
  Kokkos::View<double**, Kokkos::LayoutLeft, DeviceType> Y ("Y", numRows, numVecs);
  auto Y_h = Kokkos::create_mirror_view (Y); // we'll want this later

  // Start with the easy test case, where the matrix and the vectors
  // are all double.
  {
    constexpr double ZERO_d = static_cast<double> (0.0);
    constexpr double ONE_d = static_cast<double> (1.0);
    constexpr double TWO_d = static_cast<double> (2.0);

    double_matrix_type A_d ("A_d", G);
    auto A_d_val_h = Kokkos::create_mirror_view (A_d.values);
    A_d_val_h[0] = ONE_d;
    // This cast is deliberate; we want to use float eps here, but as
    // a double-precision number.  This is just a sanity check for
    // accuracy of the sparse mat-vec when not using mixed precision.
    A_d_val_h[1] = static_cast<double> (EPS_f) / TWO_d;
    EXPECT_NE( A_d_val_h[1], ZERO_d ); // just making sure
    Kokkos::deep_copy (A_d.values, A_d_val_h);

    // Just to make sure, we purge the previous contents of y,
    // before doing the sparse mat-vec.
    Kokkos::deep_copy (y, ZERO_d);
    KokkosSparse::spmv ("N", alpha_d, A_d, x, beta_d, y);

    Kokkos::deep_copy (y_h, y);
    const double expectedResult_allDouble = static_cast<double> (1.0) +
      static_cast<double> (EPS_f) / static_cast<double> (2.0);
    EXPECT_NE( expectedResult_allDouble, ZERO_d );
    EXPECT_EQ( y_h[0], expectedResult_allDouble );

    for (int curNumVecs = 1; curNumVecs <= numVecs; ++curNumVecs) {
      const Kokkos::pair<int, int> vecRng (0, curNumVecs);
      auto X_sub = Kokkos::subview (X, Kokkos::ALL (), vecRng);
      auto Y_sub = Kokkos::subview (Y, Kokkos::ALL (), vecRng);

      // Just to make sure, we purge the previous contents of Y,
      // before doing the sparse mat-vec.
      Kokkos::deep_copy (Y, ZERO_d);
      KokkosSparse::spmv ("N", alpha_d, A_d, X, beta_d, Y);

      Kokkos::deep_copy (Y_h, Y);
      for (int j = 0; j < curNumVecs; ++j) {
        const double actualResult = Y_h(0,j);
        EXPECT_EQ( actualResult, expectedResult_allDouble );
      }
    }
  }

  // Now exercise the case where the matrix is in float, but the
  // vectors are in double.
  {
    constexpr float ZERO_f = static_cast<float> (0.0);
    constexpr float ONE_f = static_cast<float> (1.0);
    constexpr float TWO_f = static_cast<float> (2.0);
    constexpr double ZERO_d = static_cast<double> (0.0);

    float_matrix_type A_f ("A_f", G);
    auto A_f_val_h = Kokkos::create_mirror_view (A_f.values);
    A_f_val_h[0] = ONE_f;
    A_f_val_h[1] = EPS_f / TWO_f;
    EXPECT_NE( A_f_val_h[1], ZERO_f ); // just making sure
    Kokkos::deep_copy (A_f.values, A_f_val_h);

    // Just to make sure, we purge the previous contents of y,
    // before doing the sparse mat-vec.
    Kokkos::deep_copy (y, ZERO_d);
    KokkosSparse::spmv ("N", alpha_d, A_f, x, beta_d, y);

    Kokkos::deep_copy (y_h, y);
    const double expectedResult_mixed = static_cast<double> (1.0) +
      static_cast<double> (EPS_f) / static_cast<double> (2.0);
    EXPECT_NE( expectedResult_mixed, ZERO_d );
    EXPECT_EQ( y_h[0], expectedResult_mixed );

    for (int curNumVecs = 1; curNumVecs <= numVecs; ++curNumVecs) {
      const Kokkos::pair<int, int> vecRng (0, curNumVecs);
      auto X_sub = Kokkos::subview (X, Kokkos::ALL (), vecRng);
      auto Y_sub = Kokkos::subview (Y, Kokkos::ALL (), vecRng);

      // Just to make sure, we purge the previous contents of Y,
      // before doing the sparse mat-vec.
      Kokkos::deep_copy (Y, ZERO_d);
      KokkosSparse::spmv ("N", alpha_d, A_f, X, beta_d, Y);

      Kokkos::deep_copy (Y_h, Y);
      for (int j = 0; j < curNumVecs; ++j) {
        const double actualResult = Y_h(0,j);
        EXPECT_EQ( actualResult, expectedResult_mixed );
      }
    }
  }
}

TEST( KokkosSparse, GitHubIssue101 ) {
#if defined(KOKKOS_ENABLE_CUDA)
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space> device_type;
  typedef Kokkos::View<int*, device_type>::HostMirror::device_type::execution_space host_execution_space;
  typedef Kokkos::Device<host_execution_space, Kokkos::HostSpace> host_device_type;

  test_github_issue_101<device_type> ();
  test_github_issue_101<host_device_type> ();
#else
  typedef Kokkos::View<int*>::device_type device_type;
  test_github_issue_101<device_type> ();
#endif // defined(KOKKOS_ENABLE_CUDA)
}

} // namespace (anonymous)

