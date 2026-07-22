// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Getrf.hpp>
#include <KokkosBatched_Getrs.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched getrf/getrs
/// Solving A * x = b, where
///   A: [[1, 2, 3],
///       [2, -4, 6],
///       [3, -9, -3]]
///   b: [1, 1, 1]
///   x: [23/32, 1/8, 1/96]
///   A is a real general matrix
///
/// This corresponds to the following system of equations:
///        1 x0 + 2 x1 + 3 x2 = 1
///        2 x0 - 4 x1 + 6 x2 = 1
///        3 x0 - 9 x1 - 3 x2 = 1
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType  = Kokkos::View<double**, ExecutionSpace>;
    using View3DType  = Kokkos::View<double***, ExecutionSpace>;
    using PivViewType = Kokkos::View<int**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Matrix A
    View3DType A("A", Nb, n, n);

    // Solution
    View2DType x("x", Nb, n);

    // Pivot
    PivViewType ipiv("ipiv", Nb, n);

    // Initialize A and x
    Kokkos::deep_copy(x, 1.0);
    auto h_A = Kokkos::create_mirror_view(A);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill the matrix
      h_A(ib, 0, 0) = 1.0;
      h_A(ib, 0, 1) = 2.0;
      h_A(ib, 0, 2) = 3.0;
      h_A(ib, 1, 0) = 2.0;
      h_A(ib, 1, 1) = -4.0;
      h_A(ib, 1, 2) = 6.0;
      h_A(ib, 2, 0) = 3.0;
      h_A(ib, 2, 1) = -9.0;
      h_A(ib, 2, 2) = -3.0;
    }
    Kokkos::deep_copy(A, h_A);

    // solve A * x = b with getrf and getrs
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "getrf-getrs", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A    = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_ipiv = Kokkos::subview(ipiv, ib, Kokkos::ALL);
          auto sub_x    = Kokkos::subview(x, ib, Kokkos::ALL);

          // PLU Factorize A by getrf
          KokkosBatched::SerialGetrf<KokkosBatched::Algo::Getrf::Unblocked>::invoke(sub_A, sub_ipiv);

          // Solve A * x = b with getrs
          KokkosBatched::SerialGetrs<KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Getrs::Unblocked>::invoke(
              sub_A, sub_ipiv, sub_x);
        });

    // Confirm that the results are correct
    auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - 23.0 / 32.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - 1.0 / 8.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - 1.0 / 96.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "getrf/getrs works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
