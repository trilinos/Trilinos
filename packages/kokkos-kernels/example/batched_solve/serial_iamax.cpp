// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Iamax.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched iamax
/// Finds the index of the first element having maximum absolute value.
///   X0: [1, 2, 0] -> 1
///   X1: [-5, 4, 3] -> 0
///   X2: [0, -1, -1] -> 1
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using Idx1DType  = Kokkos::View<int*, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Batched vectors
    View2DType x0("x0", Nb, n), x1("x1", Nb, n), x2("x2", Nb, n);

    // Max indices
    Idx1DType iamax0("iamax0", Nb), iamax1("iamax1", Nb), iamax2("iamax2", Nb);

    // Initialize x0, x1, x2
    auto h_x0 = Kokkos::create_mirror_view(x0);
    auto h_x1 = Kokkos::create_mirror_view(x1);
    auto h_x2 = Kokkos::create_mirror_view(x2);

    for (int ib = 0; ib < Nb; ib++) {
      h_x0(ib, 0) = 1.0;
      h_x0(ib, 1) = 2.0;
      h_x0(ib, 2) = 0.0;

      h_x1(ib, 0) = -5.0;
      h_x1(ib, 1) = 4.0;
      h_x1(ib, 2) = 3.0;

      h_x2(ib, 0) = 0.0;
      h_x2(ib, 1) = -1.0;
      h_x2(ib, 2) = -1.0;
    }
    Kokkos::deep_copy(x0, h_x0);
    Kokkos::deep_copy(x1, h_x1);
    Kokkos::deep_copy(x2, h_x2);

    // Find max indices
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "iamax", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_x0 = Kokkos::subview(x0, ib, Kokkos::ALL);
          auto sub_x1 = Kokkos::subview(x1, ib, Kokkos::ALL);
          auto sub_x2 = Kokkos::subview(x2, ib, Kokkos::ALL);

          // Find max indices
          iamax0(ib) = KokkosBatched::SerialIamax::invoke(sub_x0);
          iamax1(ib) = KokkosBatched::SerialIamax::invoke(sub_x1);
          iamax2(ib) = KokkosBatched::SerialIamax::invoke(sub_x2);
        });

    // Confirm that the results are correct
    auto h_iamax0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, iamax0);
    auto h_iamax1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, iamax1);
    auto h_iamax2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, iamax2);
    bool correct  = true;
    for (int ib = 0; ib < Nb; ib++) {
      if (h_iamax0(ib) != 1) correct = false;
      if (h_iamax1(ib) != 0) correct = false;
      if (h_iamax2(ib) != 1) correct = false;
    }

    if (correct) {
      std::cout << "iamax works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
