// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched trsm
/// Solving A * X = alpha * B, where
///   A: [[1,  1],
///       [0,  2]]
///   B: [[1, 1],
///       [1, 1]]
///   alpha: 1.5
///
///   X: [[3/4, 3/4],
///       [3/4, 3/4]]
///
/// This corresponds to the following system of equations:
///    1 x00 + 1 x10 = 1.5
///    1 x01 + 1 x11 = 1.5
///            2 x10 = 1.5
///            2 x11 = 1.5
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 2;

    // Matrix A and B
    View3DType A("A", Nb, n, n), B("B", Nb, n, n);

    // Lower triangular matrix
    // Initialize A and B
    Kokkos::deep_copy(B, 1.0);
    auto h_A = Kokkos::create_mirror_view(A);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      h_A(ib, 0, 0) = 1.0;
      h_A(ib, 0, 1) = 1.0;
      h_A(ib, 1, 0) = 0.0;
      h_A(ib, 1, 1) = 2.0;
    }
    Kokkos::deep_copy(A, h_A);

    // solve A * X = alpha * B with trsm
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "trsm", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_B = Kokkos::subview(B, ib, Kokkos::ALL, Kokkos::ALL);

          // Solve A * X = alpha * B with trsm
          KokkosBatched::SerialTrsm<KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                    KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                    KokkosBatched::Algo::Trsm::Unblocked>::invoke(alpha, sub_A, sub_B);
        });

    // Confirm that the results are correct
    auto h_B     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_B(ib, 0, 0) - 3.0 / 4.0) > eps) correct = false;
      if (Kokkos::abs(h_B(ib, 0, 1) - 3.0 / 4.0) > eps) correct = false;
      if (Kokkos::abs(h_B(ib, 1, 0) - 3.0 / 4.0) > eps) correct = false;
      if (Kokkos::abs(h_B(ib, 1, 1) - 3.0 / 4.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "trsm works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
