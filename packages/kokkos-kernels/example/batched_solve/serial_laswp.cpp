// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosBatched_Laswp.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched laswp
/// Applies a series of row interchanges to a batch of matrices.
///   A: [[4, 1, 5],
///       [2, 3, 7],
///       [6, 0, 8]]
///   p: [1, 2, 0]
///
///             Initial        0<->1          1<->2          2<->0
///   Forward:  [[4, 1, 5], -> [[2, 3, 7], -> [[2, 3, 7], -> [[4, 1, 5],
///              [2, 3, 7],     [4, 1, 5],     [6, 0, 8],     [6, 0, 8],
///              [6, 0, 8]]     [6, 0, 8]]     [4, 1, 5]]     [2, 3, 7]]
///             Initial        2<->0          1<->2          0<->1
///   Backward: [[4, 1, 5], -> [[6, 0, 8], -> [[6, 0, 8], -> [[4, 1, 5],
///              [2, 3, 7],     [2, 3, 7],     [4, 1, 5],     [6, 0, 8],
///              [6, 0, 8]]     [4, 1, 5]]     [2, 3, 7]]     [2, 3, 7]]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View3DType    = Kokkos::View<double***, ExecutionSpace>;
    using PivView2DType = Kokkos::View<int**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Batched vectors
    View3DType A("A", Nb, n, n), B("B", Nb, n, n);

    // Max indices
    PivView2DType ipiv("ipiv", Nb, n);

    // Initialize A, ipiv
    auto h_A    = Kokkos::create_mirror_view(A);
    auto h_ipiv = Kokkos::create_mirror_view(ipiv);

    for (int ib = 0; ib < Nb; ib++) {
      h_A(ib, 0, 0) = 4.0;
      h_A(ib, 0, 1) = 1.0;
      h_A(ib, 0, 2) = 5.0;
      h_A(ib, 1, 0) = 2.0;
      h_A(ib, 1, 1) = 3.0;
      h_A(ib, 1, 2) = 7.0;
      h_A(ib, 2, 0) = 6.0;
      h_A(ib, 2, 1) = 0.0;
      h_A(ib, 2, 2) = 8.0;

      h_ipiv(ib, 0) = 1;
      h_ipiv(ib, 1) = 2;
      h_ipiv(ib, 2) = 0;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(ipiv, h_ipiv);

    // Copy A into B
    Kokkos::deep_copy(B, A);

    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "laswp", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A    = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_B    = Kokkos::subview(B, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_ipiv = Kokkos::subview(ipiv, ib, Kokkos::ALL);

          // Forward direction
          KokkosBatched::SerialLaswp<KokkosBatched::Direct::Forward>::invoke(sub_ipiv, sub_A);

          // Backward direction
          KokkosBatched::SerialLaswp<KokkosBatched::Direct::Backward>::invoke(sub_ipiv, sub_B);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_A, A);
    auto h_B = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B);

    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_A(ib, 0, 0) - 4.0) > eps || Kokkos::abs(h_B(ib, 0, 0) - 4.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 0, 1) - 1.0) > eps || Kokkos::abs(h_B(ib, 0, 1) - 1.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 0, 2) - 5.0) > eps || Kokkos::abs(h_B(ib, 0, 2) - 5.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 1, 0) - 6.0) > eps || Kokkos::abs(h_B(ib, 1, 0) - 6.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 1, 1) - 0.0) > eps || Kokkos::abs(h_B(ib, 1, 1) - 0.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 1, 2) - 8.0) > eps || Kokkos::abs(h_B(ib, 1, 2) - 8.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 2, 0) - 2.0) > eps || Kokkos::abs(h_B(ib, 2, 0) - 2.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 2, 1) - 3.0) > eps || Kokkos::abs(h_B(ib, 2, 1) - 3.0) > eps) correct = false;
      if (Kokkos::abs(h_A(ib, 2, 2) - 7.0) > eps || Kokkos::abs(h_B(ib, 2, 2) - 7.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "laswp works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
