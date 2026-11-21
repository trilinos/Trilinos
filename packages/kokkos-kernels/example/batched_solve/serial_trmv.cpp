// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Trmv.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched trmv
/// x:= alpha * op(A) * x
///   A: [[1, -3, -2,  0],
///       [0,  3, -1, -2],
///       [0,  0,  9,  5],
///       [0,  0,  0, 27]]
///   x: [1, 2, 3, 4]
///   alpha: 1.5
///
/// In upper, no-transposed and non-unit case
///  x: [-16.5, -7.5, 70.5, 162]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Matrix A
    View3DType A("A", Nb, n, n);

    // Vector x
    View2DType x("x", Nb, n);

    // Upper triangular matrix
    // Initialize A and x
    auto h_A = Kokkos::create_mirror_view(A);
    auto h_x = Kokkos::create_mirror_view(x);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Set A
      h_A(ib, 0, 0) = 1.0;
      h_A(ib, 0, 1) = -3.0;
      h_A(ib, 0, 2) = -2.0;
      h_A(ib, 0, 3) = 0.0;
      h_A(ib, 1, 0) = 0.0;
      h_A(ib, 1, 1) = 3.0;
      h_A(ib, 1, 2) = -1.0;
      h_A(ib, 1, 3) = -2.0;
      h_A(ib, 2, 2) = 9.0;
      h_A(ib, 2, 3) = 5.0;
      h_A(ib, 3, 3) = 27.0;

      // Set x
      // Fill vector x
      for (int i = 0; i < n; i++) {
        h_x(ib, i) = i + 1;
      }
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(x, h_x);

    // Compute x := alpha * A * x with trmv
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "trmv", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);

          // Compute x := alpha * A * x with trmv
          KokkosBatched::SerialTrmv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                    KokkosBatched::Diag::NonUnit>::invoke(alpha, sub_A, sub_x);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_x, x);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) + 16.5) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) + 7.5) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - 70.5) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 3) - 162) > eps) correct = false;
    }

    if (correct) {
      std::cout << "trmv works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
