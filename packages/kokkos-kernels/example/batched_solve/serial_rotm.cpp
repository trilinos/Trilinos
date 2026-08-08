// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Rotm.hpp>
#include <KokkosBatched_Rotmg.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched rotm and rotmg
/// rotmg constructs the modified Givens rotation matrix H which eliminates the second component of a 2-vector [d1, d2]
/// and applies the rotation to a 2-vector [x1, y1].
/// Then apply the rotation to vectors X and Y:
/// [[X'],  = H  * [[X],
///  [Y']]          [Y]]
///
/// Usage example:
/// d1 = 1.0, d2 = 1.0, x1 = 3.0, y1 = 4.0
/// X = [3, 3, 3], Y = [4, 4, 4]
/// X' = [6.25, 6.25, 6.25], Y' = [0, 0, 0]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View1DType = Kokkos::View<double*, ExecutionSpace>;
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Values d1 and d2 we want to rotate
    View1DType d1("d1", Nb), d2("d2", Nb), x1("x1", Nb), y1("y1", Nb);

    // Vector X and Y of the rot/rotg operation
    View2DType param("param", Nb, 5);
    View2DType X("X", Nb, n), Y("Y", Nb, n);

    // Initialize d1, d2, X and Y
    const double a = 3.0, b = 4.0;
    Kokkos::deep_copy(d1, 1.0);
    Kokkos::deep_copy(d2, 1.0);
    Kokkos::deep_copy(x1, a);
    Kokkos::deep_copy(y1, b);
    Kokkos::deep_copy(X, a);
    Kokkos::deep_copy(Y, b);

    // Compute a modified Givens rotation coefficients and apply the rotation to X and Y
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "rotm/rotmg", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_d1    = Kokkos::subview(d1, ib);
          auto sub_d2    = Kokkos::subview(d2, ib);
          auto sub_x1    = Kokkos::subview(x1, ib);
          auto sub_y1    = Kokkos::subview(y1, ib);
          auto sub_param = Kokkos::subview(param, ib, Kokkos::ALL());
          KokkosBatched::Rotmg::invoke(sub_d1, sub_d2, sub_x1, sub_y1, sub_param);

          auto sub_X = Kokkos::subview(X, ib, Kokkos::ALL());
          auto sub_Y = Kokkos::subview(Y, ib, Kokkos::ALL());
          KokkosBatched::SerialRotm::invoke(sub_X, sub_Y, sub_param);
        });

    // Confirm that the results are correct
    auto h_X     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, X);
    auto h_Y     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Y);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int j = 0; j < n; j++) {
        if (Kokkos::abs(h_X(ib, j) - 6.25) > eps) correct = false;
        if (Kokkos::abs(h_Y(ib, j)) > eps) correct = false;
      }
    }

    if (correct) {
      std::cout << "rotm/rotmg works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
