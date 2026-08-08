// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Dot.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched dot product
/// computing dot = x^T * y for a batch of vectors x and y.
///
/// Usage example:
///        x: [1, 2, 3]
///        y: [4, 5, 6]
///        dot: 32
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View1DType = Kokkos::View<double*, ExecutionSpace>;
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Vector x and y
    View2DType x("x", Nb, n), y("y", Nb, n);
    View1DType dot("dot", Nb);

    // Initialize x and y
    auto h_x = Kokkos::create_mirror_view(x);
    auto h_y = Kokkos::create_mirror_view(y);
    for (int ib = 0; ib < Nb; ib++) {
      for (int j = 0; j < n; j++) {
        h_x(ib, j) = j + 1;  // x: [1, 2, 3]
        h_y(ib, j) = j + 4;  // y: [4, 5, 6]
      }
    }
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    // Compute dot = x^T * y
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "dot", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_x   = Kokkos::subview(x, ib, Kokkos::ALL());
          auto sub_y   = Kokkos::subview(y, ib, Kokkos::ALL());
          auto sub_dot = Kokkos::subview(dot, ib);
          KokkosBatched::SerialDot<KokkosBatched::Trans::Transpose, 0>::invoke(sub_x, sub_y, sub_dot);
        });

    // Confirm that the results are correct
    auto h_dot   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_dot(ib) - 32) > eps) correct = false;
    }

    if (correct) {
      std::cout << "dot works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
