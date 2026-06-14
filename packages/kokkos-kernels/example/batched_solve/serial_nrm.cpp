// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Nrm.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched nrm
/// computing nrm = ||x||_2 for a batch of vectors x.
///
/// Usage example:
///        x: [1, 3, 5]
///        nrm: sqrt(1^2 + 3^2 + 5^2) = sqrt(35)
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View1DType = Kokkos::View<double*, ExecutionSpace>;
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 3;

    // Vector x
    View2DType x("x", Nb, n);
    View1DType norm("norm", Nb);

    // Initialize x
    auto h_x = Kokkos::create_mirror_view(x);
    for (int ib = 0; ib < Nb; ib++) {
      h_x(ib, 0) = 1;
      h_x(ib, 1) = 3;
      h_x(ib, 2) = 5;
    }
    Kokkos::deep_copy(x, h_x);

    // Compute L2 norm of x
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "nrm", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_x    = Kokkos::subview(x, ib, Kokkos::ALL());
          auto sub_norm = Kokkos::subview(norm, ib);
          KokkosBatched::SerialNrm<KokkosBatched::Norm::L2>::invoke(sub_x, sub_norm);
        });

    // Confirm that the results are correct
    auto h_norm  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_norm(ib) - Kokkos::sqrt(35)) > eps) correct = false;
    }

    if (correct) {
      std::cout << "nrm works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
