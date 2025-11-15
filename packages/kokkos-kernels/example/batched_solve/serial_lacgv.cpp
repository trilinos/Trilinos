// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosBatched_Lacgv.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched lacgv
/// Conjugates a complex vector of length N
///   x: [1 + 1j, -3+2j, -2-2j,  0+1j]
///   conj(x): [1 - 1j, -3-2j, -2+2j,  0-1j]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<Kokkos::complex<double>**, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Batched vectors
    View2DType x("x", Nb, n);

    // Initialize x
    auto h_x = Kokkos::create_mirror_view(x);

    for (int ib = 0; ib < Nb; ib++) {
      h_x(ib, 0) = Kokkos::complex<double>(1.0, 1.0);
      h_x(ib, 1) = Kokkos::complex<double>(-3.0, 2.0);
      h_x(ib, 2) = Kokkos::complex<double>(-2.0, -2.0);
      h_x(ib, 3) = Kokkos::complex<double>(0.0, 1.0);
    }
    Kokkos::deep_copy(x, h_x);

    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "lacgv", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);
          KokkosBatched::SerialLacgv::invoke(sub_x);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_x, x);

    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      if (Kokkos::abs(h_x(ib, 0) - Kokkos::complex<double>(1.0, -1.0)) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 1) - Kokkos::complex<double>(-3.0, -2.0)) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 2) - Kokkos::complex<double>(-2.0, 2.0)) > eps) correct = false;
      if (Kokkos::abs(h_x(ib, 3) - Kokkos::complex<double>(0.0, -1.0)) > eps) correct = false;
    }

    if (correct) {
      std::cout << "lacgv works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
