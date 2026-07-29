// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosBatched_Rot.hpp>
#include <KokkosBatched_Rotg.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched rot and rotg
/// Compute c ans s such that
/// [[c,  s],  * [[a],  = [[r],
///  [-s, c]]     [b]]     [0]]
/// Then apply the rotation to vectors X and Y:
/// [[X'],  = [[c,  s],  * [[X],
///  [Y']]     [-s, c]]     [Y]]
///
/// Usage example:
/// a = 3, b = 4
/// X = [3, 3, 3], Y = [4, 4, 4]
/// X' = [5, 5, 5], Y' = [0, 0, 0]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View1DType = Kokkos::View<double*, ExecutionSpace>;
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Values A and B we want to rotate
    View1DType A("A", Nb), B("B", Nb), C("C", Nb), S("S", Nb);

    // Vector X and Y of the rot/rotg operation
    View2DType X("X", Nb, n), Y("Y", Nb, n);

    // Initialize A, B, X and Y
    const double a = 3.0, b = 4.0;
    Kokkos::deep_copy(A, a);
    Kokkos::deep_copy(B, b);
    Kokkos::deep_copy(X, a);
    Kokkos::deep_copy(Y, b);

    // Compute givens rotation coefficients and apply the rotation to X and Y
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "rot/rotg", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib);
          auto sub_B = Kokkos::subview(B, ib);
          auto sub_C = Kokkos::subview(C, ib);
          auto sub_S = Kokkos::subview(S, ib);
          KokkosBatched::Rotg::invoke(sub_A, sub_B, sub_C, sub_S);

          auto sub_X = Kokkos::subview(X, ib, Kokkos::ALL());
          auto sub_Y = Kokkos::subview(Y, ib, Kokkos::ALL());
          KokkosBatched::SerialRot<>::invoke(sub_X, sub_Y, sub_C(), sub_S());
        });

    // Confirm that the results are correct
    auto h_X     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, X);
    auto h_Y     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Y);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int j = 0; j < n; j++) {
        if (Kokkos::abs(h_X(ib, j) - 5) > eps) correct = false;
        if (Kokkos::abs(h_Y(ib, j)) > eps) correct = false;
      }
    }

    if (correct) {
      std::cout << "rot/rotg works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
