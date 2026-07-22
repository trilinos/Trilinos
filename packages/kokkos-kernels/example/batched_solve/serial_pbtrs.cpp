// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Pbtrf.hpp>
#include <KokkosBatched_Pbtrs.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched pbtrf/pbtrs
/// Solving A * x = b, where
///   A: [[4, 1, 0],
///       [1, 4, 1],
///       [0, 1, 4]]
///   b: [1, 1, 1]
///   x: [3/14, 1/7, 3/14]
///   A is a real symmetric (or complex Hermitian) positive definite matrix
///
/// In upper band storage,
///   Ab: [[0, 1, 1],
///        [4, 4, 4],]
/// In lower band storage,
///   Ab: [[4, 4, 4],
///        [1, 1, 0],]
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
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 3, k = 1;

    // Matrix Ab in band storage
    View3DType Ab("Ab", Nb, k + 1, n);

    // Solution
    View2DType x("x", Nb, n);

    // Initialize Ab and x
    Kokkos::deep_copy(x, 1.0);
    auto h_Ab = Kokkos::create_mirror_view(Ab);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill the diagonal elements
      h_Ab(ib, 1, 0) = 4.0;
      h_Ab(ib, 1, 1) = 4.0;
      h_Ab(ib, 1, 2) = 4.0;

      // Fill the super-diagonal elements
      h_Ab(ib, 0, 1) = 1.0;
      h_Ab(ib, 0, 2) = 1.0;
    }
    Kokkos::deep_copy(Ab, h_Ab);

    // solve A * x = b with pbtrf and pbtrs
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "pbtrf-pbtrs", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_Ab = Kokkos::subview(Ab, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x  = Kokkos::subview(x, ib, Kokkos::ALL);

          // Factorize Ab by pbtrf
          KokkosBatched::SerialPbtrf<KokkosBatched::Uplo::Upper, KokkosBatched::Algo::Pbtrf::Unblocked>::invoke(sub_Ab);

          // Solve A * x = b with pbtrs
          KokkosBatched::SerialPbtrs<KokkosBatched::Uplo::Upper, KokkosBatched::Algo::Pbtrs::Unblocked>::invoke(sub_Ab,
                                                                                                                sub_x);
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
      std::cout << "pbtrf/pbtrs works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
