// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Gbtrf.hpp>
#include <KokkosBatched_Gbtrs.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched gbtrf/gbtrs
/// Solving A * x = b, where
///   A: [[1, -3, -2,  0],
///       [-1, 1, -3, -2],
///       [2, -1,  1, -3],
///       [0,  2, -1,  1]]
///   b: [1, 1, 1, 1]
///   x: [67/81, 22/81, -40/81, -1/27]
///
/// In band storage,
///   Ab: [[0,   0,  0,  0],
///        [0,   0,  0,  0],
///        [0,   0, -2, -2],
///        [0,  -3, -3, -3],
///        [1,   1,  1,  1],
///        [-1, -1, -1,  0],
///        [2,   2,  0,  0]]
///
/// This corresponds to the following system of equations:
///      x0 - 3 x1 - 2 x2        = 1
///    - x0 +   x1 - 3 x2 - 2 x3 = 1
///    2 x0 -   x1 +   x3 - 3 x3 = 1
///           2 x1 -   x2 +   x3 = 1
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType  = Kokkos::View<double**, ExecutionSpace>;
    using View3DType  = Kokkos::View<double***, ExecutionSpace>;
    using PivViewType = Kokkos::View<int**, ExecutionSpace>;
    const int Nb = 10, n = 4, kl = 2, ku = 2;
    const int ldab = 2 * kl + ku + 1;

    // Matrix Ab in band storage
    View3DType Ab("Ab", Nb, ldab, n);

    // Solution
    View2DType x("x", Nb, n);

    // Pivot
    PivViewType ipiv("ipiv", Nb, n);

    // Initialize Ab and x
    Kokkos::deep_copy(x, 1.0);
    auto h_Ab = Kokkos::create_mirror_view(Ab);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill the band matrix
      h_Ab(ib, 2, 2) = -2;
      h_Ab(ib, 2, 3) = -2;
      h_Ab(ib, 3, 1) = -3;
      h_Ab(ib, 3, 2) = -3;
      h_Ab(ib, 3, 3) = -3;
      h_Ab(ib, 4, 0) = 1;
      h_Ab(ib, 4, 1) = 1;
      h_Ab(ib, 4, 2) = 1;
      h_Ab(ib, 4, 3) = 1;
      h_Ab(ib, 5, 0) = -1;
      h_Ab(ib, 5, 1) = -1;
      h_Ab(ib, 5, 2) = -1;
      h_Ab(ib, 6, 0) = 2;
      h_Ab(ib, 6, 1) = 2;
    }
    Kokkos::deep_copy(Ab, h_Ab);

    // solve A * x = b with gbtrf and gbtrs
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "gbtrf-gbtrs", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_Ab   = Kokkos::subview(Ab, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_ipiv = Kokkos::subview(ipiv, ib, Kokkos::ALL);
          auto sub_x    = Kokkos::subview(x, ib, Kokkos::ALL);

          // Factorize Ab by gbtrf
          KokkosBatched::SerialGbtrf<KokkosBatched::Algo::Gbtrf::Unblocked>::invoke(sub_Ab, sub_ipiv, kl, ku);

          // Solve A * x = b with gbtrs
          KokkosBatched::SerialGbtrs<KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Gbtrs::Unblocked>::invoke(
              sub_Ab, sub_ipiv, sub_x, kl, ku);
        });

    // Confirm that the results are correct
    auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - 67.0 / 81.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - 22.0 / 81.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) + 40.0 / 81.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 3) + 1.0 / 27.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "gbtrf/gbtrs works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
