// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosLapack_geqrf.hpp>
#include <KokkosLapack_gemqr.hpp>

int main(void) {
  bool correct = true;
  Kokkos::initialize();
  {
    using execution_space = Kokkos::DefaultExecutionSpace;
    using KAT             = KokkosKernels::ArithTraits<double>;

    Kokkos::View<double**, Kokkos::LayoutLeft> A("A", 3, 3), R("R", 3, 3);
    Kokkos::View<double*, Kokkos::LayoutLeft> tau("tau", 3);
    Kokkos::View<int*, Kokkos::LayoutLeft> Info("Info", 1);

    auto h_A = Kokkos::create_mirror_view(A);

    h_A(0, 0) = 12;
    h_A(0, 1) = -51;
    h_A(0, 2) = 4;
    h_A(1, 0) = 6;
    h_A(1, 1) = 167;
    h_A(1, 2) = -68;
    h_A(2, 0) = -4;
    h_A(2, 1) = 24;
    h_A(2, 2) = -41;
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(R, A);

    execution_space space{};
    KokkosLapack::geqrf(space, A, tau, Info);
    // Compute Q^T*A = Q^T*Q*R = R
    KokkosLapack::gemqr(space, "L", "T", A, tau, R, Info);
    Kokkos::fence();
    auto h_R = Kokkos::create_mirror(R);
    Kokkos::deep_copy(h_R, R);

    Kokkos::View<double**, Kokkos::DefaultHostExecutionSpace> analytic("result", 3, 3);
    analytic(0, 0) = -14;
    analytic(0, 1) = -21;
    analytic(0, 2) = 14;
    analytic(1, 0) = 0;
    analytic(1, 1) = -175;
    analytic(1, 2) = 70;
    analytic(2, 0) = 0;
    analytic(2, 1) = 0;
    analytic(2, 2) = -35;

    for (int rowIdx = 0; rowIdx < h_R.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < h_R.extent_int(1); ++colIdx) {
        if (colIdx < rowIdx) {
          if (KAT::abs(h_R(rowIdx, colIdx)) > 100 * KAT::epsilon()) {
            std::cout << "h_R(" << rowIdx << ", " << colIdx << ")=" << h_R(rowIdx, colIdx) << ", analytic is "
                      << analytic(rowIdx, colIdx) << std::endl;
            correct = false;
          }
        } else if (KAT::abs(h_R(rowIdx, colIdx) - analytic(rowIdx, colIdx)) >
                   KAT::abs(analytic(rowIdx, colIdx)) * KAT::epsilon()) {
          std::cout << "h_R(" << rowIdx << ", " << colIdx << ")=" << h_R(rowIdx, colIdx) << ", analytic is "
                    << analytic(rowIdx, colIdx) << std::endl;
          correct = false;
        }
      }
    }
  }
  Kokkos::finalize();

  if (!correct) throw std::runtime_error("KokkosLapack::geqrf is incorrect!");
  return 0;
}
