// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Swap.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched swap
/// Swapping a batch of two vectors x and y
///
/// Usage example:
///   Before: x: [1, 3, 5], y: [2, 4, 6]
///   After:  x: [2, 4, 6], y: [1, 3, 5]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Vectors x and y
    View2DType x("x", Nb, n), y("y", Nb, n);

    // Initialize x and y
    auto h_x = Kokkos::create_mirror_view(x);
    auto h_y = Kokkos::create_mirror_view(y);
    for (int ib = 0; ib < Nb; ib++) {
      h_x(ib, 0) = 1;
      h_x(ib, 1) = 3;
      h_x(ib, 2) = 5;
      h_y(ib, 0) = 2;
      h_y(ib, 1) = 4;
      h_y(ib, 2) = 6;
    }
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    // Swap x and y
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "swap", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL());
          auto sub_y = Kokkos::subview(y, ib, Kokkos::ALL());
          KokkosBatched::SerialSwap::invoke(sub_x, sub_y);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_x, x);
    Kokkos::deep_copy(h_y, y);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - 2) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - 4) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - 6) > eps) correct = false;
      if (Kokkos::abs(h_y(ib, 0) - 1) > eps) correct = false;
      if (Kokkos::abs(h_y(ib, 1) - 3) > eps) correct = false;
      if (Kokkos::abs(h_y(ib, 2) - 5) > eps) correct = false;
    }

    if (correct) {
      std::cout << "swap works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
