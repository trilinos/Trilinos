// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosLapack_geqrf.hpp>

int main(void) {
  bool correct = true;
  Kokkos::initialize();
  {
    using execution_space = Kokkos::DefaultExecutionSpace;
    using KAT             = KokkosKernels::ArithTraits<double>;

    Kokkos::View<double**, Kokkos::LayoutLeft> A("A", 3, 3);
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

    execution_space space{};
    KokkosLapack::geqrf(space, A, tau, Info);
    Kokkos::fence();
    Kokkos::deep_copy(h_A, A);

    Kokkos::View<double**, Kokkos::DefaultHostExecutionSpace> analytic("result", 3, 3);
    analytic(0, 0) = -14;
    analytic(0, 1) = -21;
    analytic(0, 2) = 14;
    analytic(1, 0) = 6. / 26;
    analytic(1, 1) = -175;
    analytic(1, 2) = 70;
    analytic(2, 0) = -4. / 26;
    analytic(2, 1) = 1. / 18;
    analytic(2, 2) = -35;

    for (int rowIdx = 0; rowIdx < A.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < A.extent_int(1); ++colIdx) {
        if (KAT::abs(h_A(rowIdx, colIdx) - analytic(rowIdx, colIdx)) >
            KAT::abs(analytic(rowIdx, colIdx)) * KAT::epsilon()) {
          std::cout << "h_A(" << rowIdx << ", " << colIdx << ")=" << h_A(rowIdx, colIdx) << ", analytic is "
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
