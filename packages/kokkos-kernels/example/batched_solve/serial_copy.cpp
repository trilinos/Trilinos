// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Copy_Decl.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched copy
/// Performs B = Op(A) where Op is one of NoTranspose, Transpose, ConjTranspose
/// A and B are 1D or 2D views
///
/// Usage example: copy A to B
/// A: [1,2,3]
/// B: [1,2,3]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Vector A and B
    View2DType A("A", Nb, n), B("B", Nb, n);

    // Initialize A and B
    auto h_A = Kokkos::create_mirror_view(A);
    auto h_B = Kokkos::create_mirror_view(B);
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vector A
      for (int i = 0; i < n; i++) {
        h_A(ib, i) = i + 1;
      }
    }
    Kokkos::deep_copy(A, h_A);

    // Compute B = A
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "copy", policy, KOKKOS_LAMBDA(int ib) {
          // B = A
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL());
          auto sub_B = Kokkos::subview(B, ib, Kokkos::ALL());
          KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose>::invoke(sub_A, sub_B);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_B, B);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_B(ib, 0) - 1) > eps) correct = false;
      if (Kokkos::abs(h_B(ib, 1) - 2) > eps) correct = false;
      if (Kokkos::abs(h_B(ib, 2) - 3) > eps) correct = false;
    }

    if (correct) {
      std::cout << "copy works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
