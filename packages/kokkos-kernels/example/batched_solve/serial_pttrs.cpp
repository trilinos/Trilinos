// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Pttrf.hpp>
#include <KokkosBatched_Pttrs.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched pttrf/pttrs
/// Solving A * x = b, where
///   A: [[4, 1, 0],
///       [1, 4, 1],
///       [0, 1, 4]]
///   b: [1, 1, 1]
///   x: [3/14, 1/7, 3/14]
/// In tridiagonal storage,
///   d: [4, 4, 4] (diagonal elements)
///   e: [1, 1] (sub/super-diagonal elements)
///
/// This corresponds to the following system of equations:
///        4 x0 +   x1        = 1
///          x0 + 4 x1 +   x2 = 1
///                 x1 + 4 x2 = 1
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Matrix A in tridiagonal storage
    View2DType d("d", Nb, n), e("e", Nb, n - 1);

    // Solution
    View2DType x("x", Nb, n);

    // Initialize d, e, and x with deep_copy
    const double d0 = 4.0, e0 = 1.0, x0 = 1.0;
    Kokkos::deep_copy(d, d0);
    Kokkos::deep_copy(e, e0);
    Kokkos::deep_copy(x, x0);

    // solve A * x = b with pttrf and pttrs
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "pttrf-pttrs", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_d = Kokkos::subview(d, ib, Kokkos::ALL);
          auto sub_e = Kokkos::subview(e, ib, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);

          // Factorize d and e by pttrf
          KokkosBatched::SerialPttrf<KokkosBatched::Algo::Pttrf::Unblocked>::invoke(sub_d, sub_e);

          // Solve A * x = b with pttrs
          KokkosBatched::SerialPttrs<KokkosBatched::Uplo::Upper, KokkosBatched::Algo::Pttrs::Unblocked>::invoke(
              sub_d, sub_e, sub_x);
        });

    // Confirm that the results are correct
    auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - 3.0 / 14.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - 1.0 / 7.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - 3.0 / 14.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "pttrf/pttrs works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
