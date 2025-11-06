// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched trsv
/// Solving A * x = alpha * b, where
///   A: [[1,  0,  0],
///       [1,  2,  0],
///       [1,  2,  3]]
///   b: [1, 1, 1]
///   alpha: 1.5
///
///   x: [3/4, 1/4, 1/2]
///
/// In lower and transposed case
///  A^T: [[1,  1,  1],
///        [0,  2,  2],
///        [0,  0,  3]]
///
/// This corresponds to the following system of equations:
///      x0 +   x1 +   x2 = 1.5
///           2 x1 + 2 x2 = 1.5
///                  3 x2 = 1.5
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Matrix A in banded storage
    View3DType A("A", Nb, n, n);

    // Vector x
    View2DType x("x", Nb, n);

    // Lower triangular matrix
    // Initialize A and x
    Kokkos::deep_copy(x, 1.0);
    auto h_A = Kokkos::create_mirror_view(A);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      h_A(ib, 0, 0) = 1.0;
      h_A(ib, 0, 1) = 0.0;
      h_A(ib, 0, 2) = 0.0;
      h_A(ib, 1, 0) = 1.0;
      h_A(ib, 1, 1) = 2.0;
      h_A(ib, 1, 2) = 0.0;
      h_A(ib, 2, 0) = 1.0;
      h_A(ib, 2, 1) = 2.0;
      h_A(ib, 2, 2) = 3.0;
    }
    Kokkos::deep_copy(A, h_A);

    // solve A^T * x = b with trsv
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "trsv", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);

          // Solve A^T * x = alpha * b with trsv
          KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose,
                                    KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(alpha,
                                                                                                                sub_A,
                                                                                                                sub_x);
        });

    // Confirm that the results are correct
    auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - 3.0 / 4.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - 1.0 / 4.0) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - 1.0 / 2.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "trsv works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
