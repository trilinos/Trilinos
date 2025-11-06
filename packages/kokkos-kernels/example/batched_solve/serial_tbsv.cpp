// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Tbsv.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched tbsv
/// Solving A * x = b, where
///   A: [[1,  1,  1],
///       [0,  2,  2],
///       [0,  0,  3]]
///   b: [1, 1, 1]
///   x: [1/2, 1/6, 1/3]
///
/// In upper band storage,
///  Ab: [[0,  0,  1],
///       [0,  1,  2],
///       [1,  2,  3]]
///
/// This corresponds to the following system of equations:
///      x0 +   x1 +   x2 = 1
///           2 x1 + 2 x2 = 1
///                  3 x2 = 1
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 3, k = 2;

    // Matrix A in banded storage
    View3DType Ab("Ab", Nb, k + 1, n);

    // Vector x
    View2DType x("x", Nb, n);

    // Upper triangular matrix
    // Initialize Ab and x
    Kokkos::deep_copy(x, 1.0);
    auto h_Ab = Kokkos::create_mirror_view(Ab);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill the diagonal elements
      h_Ab(ib, 2, 0) = 1.0;
      h_Ab(ib, 2, 1) = 2.0;
      h_Ab(ib, 2, 2) = 3.0;

      // Fill the super-diagonal elements
      h_Ab(ib, 1, 1) = 1.0;
      h_Ab(ib, 1, 2) = 2.0;

      h_Ab(ib, 0, 2) = 1.0;
    }
    Kokkos::deep_copy(Ab, h_Ab);

    // solve A * x = b with tbsv
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "tbsv", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_Ab = Kokkos::subview(Ab, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x  = Kokkos::subview(x, ib, Kokkos::ALL);

          // Solve A * x = b with tbsv
          KokkosBatched::SerialTbsv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                    KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Tbsv::Unblocked>::invoke(sub_Ab,
                                                                                                                sub_x,
                                                                                                                k);
        });

    // Confirm that the results are correct
    auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - 1.0 / 2.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - 1.0 / 6.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - 1.0 / 3.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "tbsv works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
