// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Symv.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched symv
/// Perform a symmetric matrix-vector multiplication, where
///   A: [[1,  -3, -2,  0],
///       [-3,  3, -1, -2],
///       [-2, -1,  9,  5],
///       [0,  -2,  5, 27]]
///   x: [1, 2, 3, 4], y: [5, 6, 7, 8]
///   alpha: 1.5, beta: 1.2
///
/// After, y = alpha * A * x + beta * y, it will give
///   y: [-10.5  -4.8  72.9 188.1]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Matrix A
    View3DType A("A", Nb, n, n);

    // Vector x and y
    View2DType x("x", Nb, n), y("y", Nb, n), Ref("Ref", Nb, n);

    // Initialize A, x and y
    auto h_A   = Kokkos::create_mirror_view(A);
    auto h_x   = Kokkos::create_mirror_view(x);
    auto h_y   = Kokkos::create_mirror_view(y);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vector x
      for (int i = 0; i < n; i++) {
        h_x(ib, i) = i + 1;
      }

      // Fill vector y
      for (int i = 0; i < n; i++) {
        h_y(ib, i) = 5 + i;
      }

      // Fill the matrix A
      h_A(ib, 0, 0) = 1;
      h_A(ib, 0, 1) = -3;
      h_A(ib, 0, 2) = -2;
      h_A(ib, 0, 3) = 0;
      h_A(ib, 1, 0) = -3;
      h_A(ib, 1, 1) = 3;
      h_A(ib, 1, 2) = -1;
      h_A(ib, 1, 3) = -2;
      h_A(ib, 2, 0) = -2;
      h_A(ib, 2, 1) = -1;
      h_A(ib, 2, 2) = 9;
      h_A(ib, 2, 3) = 5;
      h_A(ib, 3, 0) = 0;
      h_A(ib, 3, 1) = -2;
      h_A(ib, 3, 2) = 5;
      h_A(ib, 3, 3) = 27;

      h_Ref(ib, 0) = -10.5;
      h_Ref(ib, 1) = -4.8;
      h_Ref(ib, 2) = 72.9;
      h_Ref(ib, 3) = 188.1;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    // Compute y = alpha * A * x + beta * y
    const double alpha = 1.5, beta = 1.2;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "symv", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);
          auto sub_y = Kokkos::subview(y, ib, Kokkos::ALL);

          // y = alpha * A * x + beta * y
          KokkosBatched::SerialSymv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose>::invoke(
              alpha, sub_A, sub_x, beta, sub_y);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_y, y);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int i = 0; i < n; i++) {
        if (Kokkos::abs(h_y(ib, i) - h_Ref(ib, i)) > eps) correct = false;
      }
    }

    if (correct) {
      std::cout << "symv works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
