// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Axpy.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched axpy
/// Add entries of X scaled by coeffecient a to entries of Y: Y += a*X
///   x: [1, 2, 3, 4]
///   y: [5, 6, 7, 8]
///   a: 1.5
///
/// After, y = y + a * x, it will give
///   y: [6.5, 9.0, 11.5, 14.0]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View1DType = Kokkos::View<double*, ExecutionSpace>;
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Scalar alpha
    View1DType alpha("alpha", Nb);

    // Vector x
    View2DType x("x", Nb, n), y("y", Nb, n);

    // Initialize alpha, x and y
    Kokkos::deep_copy(alpha, 1.5);
    auto h_x = Kokkos::create_mirror_view(x);
    auto h_y = Kokkos::create_mirror_view(y);
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vector x and y
      for (int i = 0; i < n; i++) {
        h_x(ib, i) = i + 1;
        h_y(ib, i) = i + 5;
      }
    }
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    // Compute y = y + a * x
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "axpy", policy, KOKKOS_LAMBDA(int ib) {
          // y = y + a * x
          auto a     = alpha(ib);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL());
          auto sub_y = Kokkos::subview(y, ib, Kokkos::ALL());
          KokkosBatched::SerialAxpy::invoke(a, sub_x, sub_y);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_y, y);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_y(ib, 0) - 6.5) > eps) correct = false;
      if (Kokkos::abs(h_y(ib, 1) - 9.0) > eps) correct = false;
      if (Kokkos::abs(h_y(ib, 2) - 11.5) > eps) correct = false;
      if (Kokkos::abs(h_y(ib, 3) - 14.0) > eps) correct = false;
    }

    if (correct) {
      std::cout << "axpy works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
